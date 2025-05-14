import pandas as pd
import os

def filter_and_export_pfams(all_pfam_file='ALL_PFAM.txt', pfam_list_file='pfam_list.txt', output_dir='filtered_uniprot_results'):
    # 创建输出文件夹
    os.makedirs(output_dir, exist_ok=True)

    # 读取 ALL_PFAM.txt，注意编码问题
    df = pd.read_csv(all_pfam_file, sep='\t', encoding='utf-8', dtype=str)

    # 去除列名首尾空格
    df.columns = df.columns.str.strip()

    # 将需要的列转换为适当类型
    df['SSN Query Cluster #'] = pd.to_numeric(df['SSN Query Cluster #'], errors='coerce')
    df['Query-Neighbor Distance'] = pd.to_numeric(df['Query-Neighbor Distance'], errors='coerce')

    # 读取 Pfam 列表
    with open(pfam_list_file, 'r', encoding='utf-8') as f:
        pfam_list = [line.strip() for line in f if line.strip()]

    # 遍历每个 Pfam 编号
    for pfam in pfam_list:
        # 匹配 Neighbor Pfam 中包含 pfam 的行
        matched = df[df['Neighbor Pfam'].str.contains(pfam, na=False)]

        # 筛选距离小于等于10
        matched = matched[matched['Query-Neighbor Distance'] <= 10]

        # 按照 SSN Query Cluster # 升序排序
        matched = matched.sort_values(by='SSN Query Cluster #')

        # 保存到指定目录下以 pfam 命名的 CSV 文件
        if not matched.empty:
            output_path = os.path.join(output_dir, f"{pfam}.csv")
            matched.to_csv(output_path, index=False, encoding='utf-8-sig')
            print(f"[✔] {pfam} saved with {len(matched)} rows.")
        else:
            print(f"[ ] {pfam} has no matched records.")

if __name__ == "__main__":
    filter_and_export_pfams()
