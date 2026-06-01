#!/usr/bin/env python3
import os
import sys
import subprocess
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed, TimeoutError
from Bio import SeqIO
from collections import defaultdict
import shutil
import time

# --- SET PATH ---
MAFFT_BIN = "mafft"
PAL2NAL_BIN = "pal2nal.pl"
AXTCONVERT_BIN = "XTConvertor"
KAKS_BIN = "KaKs_Calculator"
KAKS_MODEL = "MA"   # Select model

def make_dict(fasta_file):
    seen = defaultdict(int)
    out = {}
    for rec in SeqIO.parse(fasta_file, "fasta"):
        seen[rec.id] += 1
        key = rec.id if seen[rec.id] == 1 else f"{rec.id}_{seen[rec.id]}"
        out[key] = rec
    return out


def run_mafft(seq1, seq2, workdir):
    prot_file = os.path.join(workdir, "prot.fa")
    SeqIO.write([seq1, seq2], prot_file, "fasta")
    aln_file = os.path.join(workdir, "prot.aln.fa")
    cmd = f"{MAFFT_BIN} --auto --thread 1 --quiet {prot_file} > {aln_file}"
    subprocess.check_call(cmd, shell=True)
    return aln_file

def run_pal2nal_clustal(prot_aln, cds1, cds2, workdir):
    nuc_file = os.path.join(workdir, "nuc.fa")
    SeqIO.write([cds1, cds2], nuc_file, "fasta")
    clustal_file = os.path.join(workdir, "codon.aln.clustal")
    cmd = f"{PAL2NAL_BIN} {prot_aln} {nuc_file} -output clustal -nogap > {clustal_file}"
    subprocess.check_call(cmd, shell=True)
    return clustal_file

def run_axtconvertor(clustal_file, workdir, tag):
    axt_file = os.path.join(workdir, f"{tag}.axt")
    cmd = [AXTCONVERT_BIN, clustal_file, axt_file]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        if os.path.exists(axt_file) and os.path.getsize(axt_file) > 0:
            pass
        else:
            return None, f"AXT Error: {proc.stderr}"
    return axt_file, None

def run_kaks(axt_file, workdir, tag):
    out_file = os.path.join(workdir, f"{tag}.kaks.txt")
    cmd = f"{KAKS_BIN} -i {axt_file} -o {out_file} -m {KAKS_MODEL}"
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if proc.returncode != 0:
        return None, f"KaKs Error: {proc.stderr}"
    return out_file, None

def parse_kaks(out_file):
    if not out_file or not os.path.exists(out_file):
        return None
    with open(out_file) as fh:
        lines = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    if len(lines) < 2:
        return None
    parts = lines[1].split()
    try:
        return {
            "seq1": parts[0],
            "seq2": parts[1],
            "Ka": parts[2],
            "Ks": parts[3],
            "Ka/Ks": parts[4],
            "P-value": parts[5],
            "S-sites": parts[6] if len(parts) > 6 else "NA",
            "N-sites": parts[7] if len(parts) > 7 else "NA"
        }
    except IndexError:
        return None

# --- worker ---
def process_pair(data):
    seq1_id, seq2_id, prot1_rec, prot2_rec, cds1_rec, cds2_rec, work_root, keep_workdirs = data
    workdir = os.path.join(work_root, f"work_{seq1_id}_{seq2_id}")
    os.makedirs(workdir, exist_ok=True)
    try:
        prot_aln = run_mafft(prot1_rec, prot2_rec, workdir)
        clustal_file = run_pal2nal_clustal(prot_aln, cds1_rec, cds2_rec, workdir)
        axt_file, err = run_axtconvertor(clustal_file, workdir, f"{seq1_id}_{seq2_id}")
        if err:
            if not keep_workdirs:
                shutil.rmtree(workdir, ignore_errors=True)
            return None, err
        kaks_out, err = run_kaks(axt_file, workdir, f"{seq1_id}_{seq2_id}")
        if err:
            if not keep_workdirs:
                shutil.rmtree(workdir, ignore_errors=True)
            return None, err
        res = parse_kaks(kaks_out)
        if res and not keep_workdirs:
            try:
                shutil.rmtree(workdir)
            except Exception:
                pass
        return res, None
    except Exception as e:
        if not keep_workdirs:
            shutil.rmtree(workdir, ignore_errors=True)
        return None, str(e)


def main():
    parser = argparse.ArgumentParser(description="Calculate Ka/Ks in parallel (bounded).")
    parser.add_argument("cds", help="CDS FASTA file")
    parser.add_argument("prot", help="Protein FASTA file")
    parser.add_argument("pairs", help="Tab-delimited pairs file (seq1 seq2)")
    parser.add_argument("output", help="Output TSV file")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of concurrent workers (default: 1)")
    parser.add_argument("-w", "--workdir", default=".", help="Root folder for per-pair work dirs")
    parser.add_argument("--keep-workdirs", action="store_true", help="Keep per-pair work directories")
    parser.add_argument("--future-timeout", type=int, default=3600, help="Seconds to wait for a worker to finish before reporting timeout")
    args = parser.parse_args()

    if not os.path.exists(AXTCONVERT_BIN):
        sys.exit(f"Error: AXTConvertor not found at {AXTCONVERT_BIN}")
    if not os.path.exists(KAKS_BIN):
        sys.exit(f"Error: KaKs_Calculator not found at {KAKS_BIN}")

    work_root = os.path.abspath(args.workdir)
    os.makedirs(work_root, exist_ok=True)

    print("Loading sequences...", flush=True)
    cds_dict = make_dict(args.cds)
    prot_dict = make_dict(args.prot)

    tasks = []
    with open(args.pairs) as pf:
        for line in pf:
            if not line.strip(): continue
            parts = line.strip().split()
            if len(parts) < 2: continue
            seq1, seq2 = parts[0], parts[1]
            if seq1 not in prot_dict or seq2 not in prot_dict:
                print(f"[Warn] {seq1} or {seq2} missing in protein FASTA. Skipping.", flush=True)
                continue
            if seq1 not in cds_dict or seq2 not in cds_dict:
                print(f"[Warn] {seq1} or {seq2} missing in CDS FASTA. Skipping.", flush=True)
                continue
            tasks.append((seq1, seq2, prot_dict[seq1], prot_dict[seq2], cds_dict[seq1], cds_dict[seq2], work_root, args.keep_workdirs))

    total = len(tasks)
    if total == 0:
        print("No tasks to run.", flush=True)
        sys.exit(0)

    print(f"Processing {total} pairs with up to {args.threads} concurrent workers (future-timeout={args.future_timeout}s)...", flush=True)

    results = []
    fails = []
    submitted = 0
    finished = 0

    # iterator over tasks
    task_iter = iter(tasks)
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        # map of future -> task ids
        future_map = {}

        # submit initial batch up to threads
        for _ in range(min(args.threads, total)):
            try:
                task = next(task_iter)
            except StopIteration:
                break
            future = executor.submit(process_pair, task)
            future_map[future] = (task[0], task[1])
            submitted += 1

        try:
            while future_map:
                # wait for any future to complete with a reasonable timeout loop
                done, _ = as_completed(future_map, timeout=args.future_timeout).__next__(), None
                # above line yields one completed future; handle it
                future = done
                s1, s2 = future_map.pop(future)
                finished += 1
                try:
                    res, err = future.result(timeout=1)
                except TimeoutError:
                    # the future is completed per as_completed but result timed out briefly; try again
                    try:
                        res, err = future.result()
                    except Exception as e:
                        res = None
                        err = f"Error retrieving result: {e}"
                except Exception as e:
                    res = None
                    err = str(e)

                if res:
                    results.append(res)
                    print(f"[{finished}/{total}] OK {s1} vs {s2} Ka/Ks={res.get('Ka/Ks','NA')}", flush=True)
                else:
                    fails.append((s1, s2, err if err else "unknown"))
                    print(f"[{finished}/{total}] FAIL {s1} vs {s2}: {err}", flush=True)

                # submit next task if available
                try:
                    task = next(task_iter)
                    future = executor.submit(process_pair, task)
                    future_map[future] = (task[0], task[1])
                    submitted += 1
                except StopIteration:
                    pass

        except StopIteration:
            # 
            pass
        except TimeoutError:
            # 
            print("Timeout: no tasks finished within the specified future-timeout.", flush=True)
            # 
            for fut in future_map:
                fut.cancel()
            print("Cancelled remaining tasks. Exiting.", flush=True)
            executor.shutdown(wait=False)
            sys.exit(1)
        except KeyboardInterrupt:
            print("Interrupted by user; cancelling running tasks...", flush=True)
            for fut in future_map:
                fut.cancel()
            executor.shutdown(wait=False)
            sys.exit(1)

    # write outputs
    if results:
        with open(args.output, "w") as out:
            header = ["seq1","seq2","Ka","Ks","Ka/Ks","P-value","S-sites","N-sites"]
            out.write("\t".join(header) + "\n")
            for r in results:
                out.write("\t".join(r[h] for h in header) + "\n")
        print(f"Written {len(results)} results to {args.output}", flush=True)
    else:
        print("No successful results to write.", flush=True)

    if fails:
        failed_file = os.path.splitext(args.output)[0] + ".failed.tsv"
        with open(failed_file, "w") as ff:
            ff.write("seq1\tseq2\treason\n")
            for s1, s2, reason in fails:
                ff.write(f"{s1}\t{s2}\t{reason}\n")
        print(f"Logged {len(fails)} failed pairs to {failed_file}", flush=True)

if __name__ == "__main__":
    main()
