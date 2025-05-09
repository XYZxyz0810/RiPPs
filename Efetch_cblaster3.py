#!/usr/bin/python
# -*- coding: UTF-8 -*-

import time
import logging
import os
from Bio import Entrez
from Bio import SeqIO
import re
from concurrent.futures import ThreadPoolExecutor, TimeoutError
from typing import List

# 设置邮箱
Entrez.email = "1157242839@qq.com"

# 配置日志记录
logging.basicConfig(filename='fetch_data.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def check_existing_record(nucleotide_file, Record_id):
    if not os.path.exists(nucleotide_file):
        return False
    with open(nucleotide_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(f">{Record_id}"):
                return True
    return False


def read_accession_list(accession_name):
    with open(accession_name, "r") as file:
        accession_lines = [line.strip() for line in file if line.strip()]
    logging.info(f"Read {len(accession_lines)} lines from accession file {accession_name}.")
    return accession_lines


def efetch_request(ids: List[str]):
    try:
        return Entrez.efetch(
            db="ipg",
            rettype="ipg",
            retmode="text",
            id=ids,
            retmax=10000,
        )
    except Exception as e:
        logging.error("Network error while retrieving IPG info: %s", e)
        raise


def efetch_IPGs(ids):
    table = []
    for start in range(0, len(ids), 10000):
        chunk = ids[start: start + 10000]
        handle = efetch_request(chunk)
        if handle.code != 200:
            raise RuntimeError(f"Bad response from NCBI [code {handle.code}]")
        for line in handle:
            if isinstance(line, bytes):
                line = line.decode()
            line = line.strip('\n')
            table.append(line)
    return table


def parse_ipg_table(ipg_lines):
    ipg_dict = {}
    header = ipg_lines[0].split('\t')
    for row in ipg_lines[1:]:
        parts = row.split('\t')
        if len(parts) < 11:
            continue
        protein = parts[6]
        nucleotide = parts[2]
        start = parts[3]
        stop = parts[4]
        strand = parts[5]

        if protein in ipg_dict:
            continue
        try:
            start_int = int(start)
            stop_int = int(stop)
        except ValueError:
            continue

        if nucleotide:
            ipg_dict[protein] = {
                "nucleotide": nucleotide,
                "start": start_int,
                "stop": stop_int,
                "strand": strand
            }
    return ipg_dict


def parse_cblaster_table(input_file):
    cblaster_dic = {}
    with open(input_file, "r", encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            parts = re.split(r'\s{2,}|\t+', line.strip())
            if len(parts) < 5:
                continue
            organism = parts[0]
            nucleotide = parts[1]
            try:
                start = int(parts[2])
                stop = int(parts[3])
            except ValueError:
                continue

            cblaster_dic[nucleotide] = {
                "nucleotide": nucleotide,
                "start": start,
                "stop": stop,
            }
    logging.info(f"Parsed {len(cblaster_dic)} records from cblaster table {input_file}.")
    return cblaster_dic

def fetch_nucleotide_sequence(genome_id, start, stop):
    try:
        handle = Entrez.efetch(db="nuccore", rettype="gb", id=genome_id, seq_start=start, seq_stop=stop)
        record = SeqIO.read(handle, format="gb")
        handle.close()
        return record
    except Exception as e:
        logging.warning(f"Fetch nucleotide sequence error: {e}")
        time.sleep(5)
        return fetch_nucleotide_sequence(genome_id, start, stop)


def process_accession_line(idx, Record_id, UDlength, nucleotide_file, failfilename, accession_lines, ipg_dict):
    logging.info(f"Processing {idx}/{len(accession_lines)}...{Record_id}...")

    if check_existing_record(nucleotide_file, Record_id):
        logging.info(f"Record {Record_id} already exists, skipping.")
        return

    if Record_id not in ipg_dict:
        with open(failfilename, "a") as f:
            f.write(f"{Record_id} not found in IPG\n")
        logging.warning(f"{Record_id} not found in IPG data.")
        return

    info = ipg_dict[Record_id]
    GenomeGI = info["nucleotide"]
    Gene_up = info["start"]
    Gene_down = info["stop"]
    strand = info["strand"]

    Nucl_up = max(0, Gene_up - UDlength) if Gene_up > UDlength else 0
    Nucl_down = Gene_down + UDlength

    NuclRecord = fetch_nucleotide_sequence(GenomeGI, Nucl_up, Nucl_down)

    if strand == "+":
        seq = str(NuclRecord.seq)
    else:
        seq = str(NuclRecord.seq.reverse_complement())

    with open(nucleotide_file, "a") as f:
        f.write(f">{Record_id} {Nucl_up}...{Nucl_down}\n{seq}\n")

    logging.info(f"Saved {Record_id}: {GenomeGI} ({Nucl_up}...{Nucl_down})")
    time.sleep(0.5)


def process_cblaster_line(idx, nucleotide_id, UDlength, nucleotide_file, failfilename, accession_lines, cblaster_dict):
    logging.info(f"Processing {idx}/{len(accession_lines)}...{nucleotide_id}...")

    if check_existing_record(nucleotide_file, nucleotide_id):
        logging.info(f"Record {nucleotide_id} already exists, skipping.")
        return

    if nucleotide_id not in cblaster_dict:
        with open(failfilename, "a") as f:
            f.write(f"{nucleotide_id} not found in cblaster table\n")
        logging.warning(f"{nucleotide_id} not found in cblaster data.")
        return

    info = cblaster_dict[nucleotide_id]
    GenomeGI = info["nucleotide"]
    Gene_up = info["start"]
    Gene_down = info["stop"]

    Nucl_up = max(0, Gene_up - UDlength) if Gene_up > UDlength else 0
    Nucl_down = Gene_down + UDlength

    NuclRecord = fetch_nucleotide_sequence(GenomeGI, Nucl_up, Nucl_down)

    seq = str(NuclRecord.seq)  # 默认正链

    with open(nucleotide_file, "a") as f:
        f.write(f">{nucleotide_id} {Nucl_up}...{Nucl_down}\n{seq}\n")

    logging.info(f"Saved {nucleotide_id}: {GenomeGI} ({Nucl_up}...{Nucl_down})")
    time.sleep(0.5)


def run_fetch_nucl(accession_name, UDlength):
    UDstr = f"{int(UDlength // 1000)}kb" if UDlength > 1000 else f"{UDlength}bp"
    Nuclfile = f"{accession_name.split('.')[0]}_{UDstr}_nucl.fasta"
    failfilename = f"{accession_name.split('.')[0]}_fail_list.txt"

    accession_lines = read_accession_list(accession_name)
    id_list = accession_lines
    ipg_table = efetch_IPGs(id_list)

    ipg_output_file = f"{accession_name.split('.')[0]}_ipg_table.txt"
    with open(ipg_output_file, "w", encoding="utf-8") as f:
        f.write("\n".join(ipg_table))
    logging.info(f"IPG table saved to {ipg_output_file}")

    ipg_dict = parse_ipg_table(ipg_table)

    ipg_txt_file = f"{accession_name.split('.')[0]}_ipg_dict.txt"
    with open(ipg_txt_file, "w", encoding="utf-8") as f:
        for protein, info in ipg_dict.items():
            line = f"{protein}\t{info['nucleotide']}\t{info['start']}\t{info['stop']}\t{info['strand']}\n"
            f.write(line)

    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {
            executor.submit(process_accession_line, idx, record_id, UDlength, Nuclfile, failfilename, accession_lines, ipg_dict): idx
            for idx, record_id in enumerate(accession_lines, 1)
        }

        for future in futures:
            try:
                future.result(timeout=30)
            except TimeoutError:
                logging.error(f"任务 {futures[future]} 超时")
            except Exception as e:
                logging.error(f"任务 {futures[future]} 异常: {e}")

    logging.info(f"Completed processing {accession_name}. Output: {Nuclfile}")
    return Nuclfile

def run_fetch_nucl_cblaster(cblaster_table, UDlength):
    UDstr = f"{int(UDlength // 1000)}kb" if UDlength > 1000 else f"{UDlength}bp"
    Nuclfile = f"{os.path.splitext(cblaster_table)[0]}_{UDstr}_nucl.fasta"
    failfilename = f"{os.path.splitext(cblaster_table)[0]}_fail_list.txt"

    cblaster_dict = parse_cblaster_table(cblaster_table)

    cblaster_output_file = f"{os.path.splitext(cblaster_table)[0]}_cblaster_dict.txt"
    with open(cblaster_output_file, "w", encoding="utf-8") as f:
        for nid, info in cblaster_dict.items():
            line = f"{nid}\t{info['nucleotide']}\t{info['start']}\t{info['stop']}\n"
            f.write(line)
    logging.info(f"cblaster dict saved to {cblaster_output_file}")

    accession_lines = list(cblaster_dict.keys())

    logging.info(f"Starting fetch for {len(accession_lines)} cblaster records...")

    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {
            executor.submit(
                process_cblaster_line,
                idx, record_id, UDlength,
                Nuclfile, failfilename,
                accession_lines, cblaster_dict
            ): idx
            for idx, record_id in enumerate(accession_lines, 1)
        }

        for future in futures:
            try:
                future.result(timeout=30)
            except TimeoutError:
                logging.error(f"任务 {futures[future]} 超时")
            except Exception as e:
                logging.error(f"任务 {futures[future]} 异常: {e}")

    logging.info(f"Completed processing {cblaster_table}. Output: {Nuclfile}")
    return Nuclfile

if __name__ == "__main__":
    accession_name = input("请输入 accession 文件名: ").strip()
    UDlength = int(input("请输入上游和下游 DNA 长度: ").strip())
    run_fetch_nucl(accession_name, UDlength)
