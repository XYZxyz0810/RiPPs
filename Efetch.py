#!/usr/bin/python
# -*- coding: UTF-8 -*-

import time
import logging
import os
from Bio import Entrez
from Bio import SeqIO
import re
import requests
from concurrent.futures import ThreadPoolExecutor, TimeoutError

# 设置邮箱
Entrez.email = "1157242839@qq.com"

# 配置日志记录
logging.basicConfig(filename='fetch_data.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# 检查输出文件是否已经包含相应的 accession 记录
def check_existing_record(nucleotide_file, Uniprot_id, Record_id):
    with open(nucleotide_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(f">{Uniprot_id}|{Record_id}"):
                return True  # 如果发现已存在记录
    return False

# 读取 accession 列表文件
def read_accession_list(accession_name):
    with open(accession_name, "r") as file:
        accession_lines = file.readlines()
    logging.info(f"Read {len(accession_lines)} lines from accession file {accession_name}.")
    return accession_lines

# 重试机制获取蛋白质记录
def fetch_with_retry(record_id, retries=3, delay=5):
    attempt = 0
    while attempt < retries:
        try:
            Prot_handle = Entrez.efetch(db="protein", rettype="gb", id=record_id)
            ProtRecord = SeqIO.read(Prot_handle, format='gb')
            Prot_handle.close()
            return ProtRecord
        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 429:
                logging.warning(f"请求过于频繁，等待 {delay} 秒后重试...")
                time.sleep(delay)
                attempt += 1
            else:
                raise
    logging.error(f"多次尝试后仍无法获取数据，跳过 ID {record_id}")
    return None

# 获取基因组信息
def get_genome_info(ProtRecord, UDlength):
    Genomelocation = str(ProtRecord.features[-1].qualifiers["coded_by"])
    NuclBar = Genomelocation.split("'")[1]
    GenomeGI = NuclBar.split("complement(")[-1].split("(")[-1].split(".")[0]
    Gene_up = int(re.findall(r"\d+", NuclBar.split(":")[1].split(".")[0])[0])
    Gene_down = int(re.findall(r"\d+", NuclBar.split(")")[0].split("..")[1])[0])
    Nucl_up = max(0, Gene_up - UDlength) if Gene_up > UDlength else 0
    Nucl_down = Gene_down + UDlength
    return GenomeGI, Nucl_up, Nucl_down

# 获取核酸序列
def fetch_nucleotide_sequence(GenomeGI, Nucl_up, Nucl_down):
    try:
        Nucl_handle = Entrez.efetch(db="nuccore", rettype="gb", id=GenomeGI, retmax="100000", seq_start=Nucl_up,
                                    seq_stop=Nucl_down)
        NuclRecord = SeqIO.read(Nucl_handle, format='gb')
        Nucl_handle.close()
        return NuclRecord
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 429:
            logging.warning(f"请求过于频繁，等待 5 秒后重试...")
            time.sleep(5)
            return fetch_nucleotide_sequence(GenomeGI, Nucl_up, Nucl_down)
        else:
            raise

# 处理 accession 行
def process_accession_line(idx, accession_line, UDlength, nucleotide_file, failfilename, accession_lines):
    Uniprot_id, Record_id = accession_line.strip().split("\t")[:2]
    logging.info(f"Processing {idx}/{len(accession_lines) - 1}...{Uniprot_id}@{Record_id}...")

    # 如果输出文件中已存在该记录，则跳过
    if check_existing_record(nucleotide_file, Uniprot_id, Record_id):
        logging.info(f"Record {Uniprot_id}@{Record_id} already exists in output file, skipping.")
        return

    ProtRecord = fetch_with_retry(Record_id)
    if ProtRecord is None:
        with open(failfilename, "a") as f:
            f.write(f"{Uniprot_id} // {Record_id} failed to fetch.\n")
        logging.warning(f"Failed to fetch {Uniprot_id}@{Record_id}, skipping.")
        return

    if ProtRecord.features[-1].type != "CDS":
        with open(failfilename, "a") as f:
            f.write(f"{Uniprot_id} // {Record_id} has no CDS.\n")
        logging.warning(f"No CDS found for {Uniprot_id}@{Record_id}, skipping.")
        return

    GenomeGI, Nucl_up, Nucl_down = get_genome_info(ProtRecord, UDlength)
    NuclRecord = fetch_nucleotide_sequence(GenomeGI, Nucl_up, Nucl_down)

    if ProtRecord.features[-1].location.strand == 1:
        seq = str(NuclRecord.seq)
    else:
        seq = str(NuclRecord.seq.reverse_complement())

    with open(nucleotide_file, "a") as f:
        f.write(f">{Uniprot_id}|{Record_id} {Nucl_up}...{Nucl_down}\n{seq}\n")

    logging.info(f"Saved sequence for {Uniprot_id}@{Record_id}: {GenomeGI} ({Nucl_up}...{Nucl_down})")
    time.sleep(0.5)

# 主处理函数
def run_fetch_nucl(accession_name, UDlength):
    UDstr = f"{int(UDlength // 1000)}kb" if UDlength > 1000 else f"{UDlength}bp"
    Nuclfile = f"{accession_name.split('.')[0]}_{UDstr}_nucl.fasta"
    failfilename = f"{accession_name.split('.')[0]}_fail_list.txt"

    accession_lines = read_accession_list(accession_name)

    with open(Nuclfile, "a") as nucleotide_file:
        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = {
                executor.submit(process_accession_line, idx, line, UDlength, Nuclfile, failfilename, accession_lines): idx
                for idx, line in enumerate(accession_lines[1:], 1)
            }

        for future in futures:
            try:
                future.result(timeout=30)
            except TimeoutError:
                logging.error(f"任务 {futures[future]} 超时")

    logging.info(f"Completed processing {accession_name}. Output: {Nuclfile}")
    return Nuclfile

# 主函数调用
if __name__ == "__main__":
    accession_name = input("请输入 accession 文件名: ").strip()
    UDlength = int(input("请输入上游和下游 DNA 长度: ").strip())
    run_fetch_data(accession_name, UDlength)
