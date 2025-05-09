import os
import time
import subprocess
import logging
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker as lm
import concurrent.futures
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from networkx.algorithms.components.connected import connected_components
from concurrent.futures import ThreadPoolExecutor  # 替换 ProcessPoolExecutor

def run_mmseqs2(peptide_file, out_dir):
    """
    运行 MMseqs2 进行蛋白质序列聚类并将结果保存至指定目录

    :param peptide_file: 输入的 peptide 文件路径
    :param out_dir: 输出目录
    :return: 聚类结果文件的路径
    """
    # 确保输出文件夹 cluster 存在
    cluster_dir = os.path.join(out_dir, "cluster")
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)
    cluster_res = os.path.join(cluster_dir, "clusterRes")
    tmp_dir = os.path.join(cluster_dir, "tmp")
    # 确保临时文件夹 tmp 存在
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    # 构建 mmseqs2 命令
    mmseqs_cmd = [
        "mmseqs", "easy-cluster", peptide_file, cluster_res, tmp_dir,
        "--min-seq-id", "0.3", "-c", "0.6", "--cov-mode", "1"
    ]
    logging.info(f"Running MMseqs2 on {peptide_file}...")
    subprocess.run(mmseqs_cmd, check=True)
    logging.info(f"MMseqs2 clustering results saved to {cluster_res}")
    # 定义聚类结果文件路径
    cluster_result_file = os.path.join(cluster_dir, "clusterRes_cluster.tsv")
    # 返回聚类结果文件路径
    return cluster_result_file


def parse_mmseqs_clusters(cluster_result_file, peptide_output_file, out_dir, min_nodes=3):
    """
    解析 MMseqs2 生成的聚类结果 TSV 文件，并生成对应的 FASTA 文件。

    :param cluster_result_file: MMseqs2 生成的聚类结果 TSV 文件
    :param peptide_file: 原始 FASTA 文件
    :param out_dir: 输出目录
    :param min_nodes: 最小的簇大小（默认 1）
    :return: 生成的 FASTA 文件路径列表
    """
    # 创建输出文件夹 "cluster_parser" 目录
    cluster_parser_dir = os.path.join(out_dir, "cluster", "cluster_parser")
    if not os.path.exists(cluster_parser_dir):
        os.makedirs(cluster_parser_dir)

    def gen(cluster_list):
        with open(cluster_list, 'r') as f:
            cont = [line for line in f]
        s = ''.join(cont).replace('\n\t', '\t').split('\n')
        return [sub.split('\t') for sub in s if len(sub) != 0]

    def to_graph(list_of_connection):
        G = nx.Graph()
        for part in list_of_connection:
            G.add_nodes_from(part)
            G.add_edges_from(to_edges(part))
        return G

    def to_edges(list_of_connection):
        it = iter(list_of_connection)
        last = next(it)
        for current in it:
            yield last, current
            last = current

    # 记录程序开始时间
    start = time.perf_counter()

    # 解析 TSV 文件
    pair = gen(cluster_result_file)
    G = to_graph(pair)
    cluster_list = list(connected_components(G))
    sorted_cluster = sorted(cluster_list, key=lambda a: len(a), reverse=True)

    # 读取 FASTA 文件
    id2seq = {seq_record.id: str(seq_record.seq) for seq_record in SeqIO.parse(peptide_output_file, 'fasta')}

    # 创建存储生成的 FASTA 文件路径列表
    n = 0
    min_nodes -= 1  # 由于索引从 0 开始，调整最小节点数
    fasta_files = []

    # 遍历聚类结果并生成对应的 FASTA 文件
    for cluster in sorted_cluster:
        if len(cluster) > min_nodes:
            n += 1
            output_fasta = os.path.join(cluster_parser_dir,
                                        f"{str(n).zfill(5)}_{len(cluster)}_{list(cluster)[0]}.fasta")
            fasta_files.append(output_fasta)
            with open(output_fasta, 'w') as f_out:
                for seq_id in cluster:
                    f_out.write(f'>{seq_id}\n{id2seq[seq_id]}\n')

    # 记录程序结束时间
    finish = time.perf_counter()
    print(f'Clustering results parsed in {round(finish - start, 3)} seconds')

    # 返回生成的 FASTA 文件路径列表
    return fasta_files


def makelogo(fasta_files, out_dir):
    # 在 cluster 目录下创建 logo 文件夹
    cluster_dir = os.path.join(out_dir, "cluster")
    logo_dir = os.path.join(cluster_dir, "logo")
    if not os.path.exists(logo_dir):
        os.makedirs(logo_dir)

    # 检查 fasta_files 是否是列表类型
    if isinstance(fasta_files, list):
        # 如果是列表，迭代每个文件处理
        for fasta_file in fasta_files:
            # 处理文件名，避免 `|` 影响 mafft
            safe_fasta = os.path.join(logo_dir, os.path.basename(fasta_file).replace("|", "_"))
            os.rename(fasta_file, safe_fasta)  # 先重命名，避免 mafft 解析错误

            # 使用 MAFFT 进行序列对齐
            mafft_cline = MafftCommandline(input=safe_fasta)
            stdout, stderr = mafft_cline()

            aligned_path = os.path.join(logo_dir, os.path.basename(safe_fasta) + ".aln")
            with open(aligned_path, "w") as handle:
                handle.write(stdout)

            # 读取对齐后的序列
            align = AlignIO.read(aligned_path, "fasta")
            seqs = [str(seq_record.seq) for seq_record in align]

            # 将对齐的序列转换为矩阵，并生成 Logo 图
            pre_aln_df = lm.alignment_to_matrix(sequences=seqs, to_type='counts', characters_to_ignore='.-X')
            pre_aln_logo = lm.Logo(pre_aln_df, font_name='sans',
                                   fade_below=0.5, shade_below=0.5,
                                   figsize=(25, 2),
                                   stack_order='big_on_top',
                                   color_scheme='NajafabadiEtAl2017',
                                   fade_probabilities=True,
                                   baseline_width=1.5,
                                   show_spines=False)

            # 保存 logo 图像
            plt.savefig(os.path.join(logo_dir, os.path.basename(safe_fasta) + ".png"))
            plt.savefig(os.path.join(logo_dir, os.path.basename(safe_fasta) + ".svg"))
            plt.close()

            # 还原原始文件名
            os.rename(safe_fasta, fasta_file)

    else:
        # 如果只有一个文件，按原来的方式处理
        safe_fasta = os.path.join(logo_dir, os.path.basename(fasta_files).replace("|", "_"))
        os.rename(fasta_files, safe_fasta)  # 先重命名，避免 mafft 解析错误

        # 使用 MAFFT 进行序列对齐
        mafft_cline = MafftCommandline(input=safe_fasta)
        stdout, stderr = mafft_cline()

        aligned_path = os.path.join(logo_dir, os.path.basename(safe_fasta) + ".aln")
        with open(aligned_path, "w") as handle:
            handle.write(stdout)

        # 读取对齐后的序列
        align = AlignIO.read(aligned_path, "fasta")
        seqs = [str(seq_record.seq) for seq_record in align]

        # 将对齐的序列转换为矩阵，并生成 Logo 图
        pre_aln_df = lm.alignment_to_matrix(sequences=seqs, to_type='counts', characters_to_ignore='.-X')
        pre_aln_logo = lm.Logo(pre_aln_df, font_name='sans',
                               fade_below=0.5, shade_below=0.5,
                               figsize=(25, 2),
                               stack_order='big_on_top',
                               color_scheme='NajafabadiEtAl2017',
                               fade_probabilities=False,
                               baseline_width=1.5,
                               show_spines=False)

        # 保存 logo 图像
        plt.savefig(os.path.join(logo_dir, os.path.basename(safe_fasta) + ".png"))
        plt.savefig(os.path.join(logo_dir, os.path.basename(safe_fasta) + ".svg"))
        plt.close()

        # 还原原始文件名
        os.rename(safe_fasta, fasta_files)



if __name__ == "__main__":
    start = time.perf_counter()  # 添加开始时间记录
    peptide_file = "Peptide_file.txt"
    output_dir = "output_dir"

    cluster_result = run_mmseqs2(peptide_file, output_dir)
    tsv_file = os.path.join(output_dir, "cluster", "clusterRes_cluster.tsv")

    fasta_files = parse_mmseqs_clusters(tsv_file, peptide_file, output_dir, min_nodes=5)

    cluster_parser_dir = os.path.join(output_dir, "cluster", "cluster_parser")
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for fasta_file in fasta_files:
            print(f'Processing {fasta_file}')
            result = executor.submit(makelogo, fasta_file, cluster_parser_dir)
            result.result()  # 等待任务完成
            print(f'Finished {fasta_file}')

    finish = time.perf_counter()  # 计算完成时间
    print(f'Finished all processes in {round(finish - start, 3)} seconds')
