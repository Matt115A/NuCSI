# run_fastp.py
#!/usr/bin/env python3
import os
import re
import glob
import yaml
import logging
import subprocess
import argparse
from datetime import datetime
from tqdm import tqdm


def load_config(path):
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    input_dir  = cfg["input_dir_raw"]
    output_dir = cfg["output_dir_qc_reads"]
    q_thr      = cfg.get("quality_threshold", 20)
    log_base   = cfg.get("log_dir_qc", os.path.join(output_dir, "logs"))
    return input_dir, output_dir, q_thr, log_base


def setup_logging(log_base):
    now     = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = os.path.join(log_base, now)
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, "run_fastp.log.txt")

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s %(levelname)s: %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    logger.info(f"Logging → {log_path}")
    return log_dir, log_path


def prepend_log(log_path, summary):
    with open(log_path) as fh:
        old = fh.read()
    with open(log_path, "w") as fh:
        fh.write(summary + "\n\n" + old)


def find_pairs(input_dir):
    pattern = os.path.join(input_dir, "*R1*.fastq.gz")
    r1_files = glob.glob(pattern)
    pairs = []
    for r1 in r1_files:
        r2 = r1.replace("R1", "R2", 1)
        if os.path.exists(r2):
            base = re.sub(r"_R1.*\.fastq\.gz$", "", os.path.basename(r1))
            pairs.append((r1, r2, base))
        else:
            logging.warning(f"Skipping {r1!r}: no matching R2 found.")
    return pairs


def run_fastp(pairs, out_dir, log_dir, q_thr, totals):
    os.makedirs(out_dir, exist_ok=True)
    for r1, r2, base in tqdm(pairs, desc="Running fastp QC"):
        out_r1 = os.path.join(out_dir, f"{base}_R1_qc.fastq.gz")
        out_r2 = os.path.join(out_dir, f"{base}_R2_qc.fastq.gz")
        json_report = os.path.join(log_dir, f"{base}_fastp.json")
        html_report = os.path.join(log_dir, f"{base}_fastp.html")

        cmd = [
            "fastp",
            "-i", r1,
            "-I", r2,
            "-o", out_r1,
            "-O", out_r2,
            "-q", str(q_thr),
            "--json", json_report,
            "--html", html_report,
            "--thread", "4",
            "--report_title", base
        ]
        logging.info(f"> {' '.join(cmd)}")
        try:
            res = subprocess.run(
                cmd, check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            logging.debug(res.stdout)
            logging.debug(res.stderr)
            # Parse the total reads from the output
            # Look for "total reads: X" pattern in the output
            m = re.search(r"total reads:\s+(\d+)", res.stdout)
            if not m:
                m = re.search(r"total reads:\s+(\d+)", res.stderr)
            got = int(m.group(1)) if m else 0
            totals["processed"] += got
            logging.info(f"[OK] {base}: {got} reads processed → {out_r1}, {out_r2}")
        except subprocess.CalledProcessError as e:
            logging.error(f"[FAIL] {base} (code {e.returncode})")
            logging.error(e.stderr.strip())


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-c","--config", required=True, help="Path to configs.yaml")
    args = p.parse_args()

    input_dir, out_dir, q_thr, log_base = load_config(args.config)
    log_dir, log_path = setup_logging(log_base)

    logging.info(f"Raw reads in:         {input_dir}")
    logging.info(f"QC outputs to:        {out_dir}")
    logging.info(f"Quality threshold:    {q_thr}")

    pairs = find_pairs(input_dir)
    logging.info(f"Found {len(pairs)} R1/R2 pairs")
    if not pairs:
        logging.warning("No pairs to process—exiting.")
        return

    totals = {"processed": 0}
    run_fastp(pairs, out_dir, log_dir, q_thr, totals)

    summary = f"TL;DR — fastp processed {totals['processed']} reads across {len(pairs)} pairs"
    logging.info(summary)
    prepend_log(log_path, summary)

if __name__ == "__main__":
    main()