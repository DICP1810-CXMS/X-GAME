import pandas as pd
import os
import re
import numpy as np


# 文件夹路径
folder_path = "/mnt/Stringent filtering of cross-links/data/TDS"

# 用于存储每个文件的 data_peptide、data_site
all_data_peptide_intra = []
all_data_site_intra = []
all_data_peptide_inter = []
all_data_site_inter = []
# 用于存储所有文件名
all_csv_filenames = []

# 遍历文件夹下的所有文件
for filename in os.listdir(folder_path):
    if filename.endswith(".csv"):
        # 构建文件路径
        file_path = os.path.join(folder_path, filename)
        all_csv_filenames.append(filename)
        
        # 读取数据
        data = pd.read_csv(file_path)
        data = data.iloc[:, 1:]  # 去掉第一列
        data_intra = data[(data['Protein_Type'] == 'Intra-Protein') | (data['Protein_Type'] == 'None')]
        data_inter = data[data['Protein_Type'] == 'Inter-Protein']
        
        if not data_intra.empty:
            PPI_library_intra = data_intra['Proteins'].str.split('/', expand=True).stack().str.strip().reset_index(drop=True)
            PPI_library_intra = PPI_library_intra[PPI_library_intra != '']
        if not data_inter.empty:
            PPI_library_inter = data_inter['Proteins'].str.split('/', expand=True).stack().str.strip().reset_index(drop=True)
            PPI_library_inter = PPI_library_inter[PPI_library_inter != '']
        
        # 定义正则表达式模式
        pattern = r'sp\|(.+)\|.+\((\d+)\)-sp\|(.+)\|.+\((\d+)\)'
        
        # 处理 intra 数据
        if not data_intra.empty:
            PPI_library_intra_filter = PPI_library_intra.str.extract(pattern)
            PPI_library_intra_filter.columns = ['Protein1', 'Site1', 'Protein2', 'Site2']
            PPI_library_intra_filter['Pair'] = PPI_library_intra_filter['Protein1'] + '-' + PPI_library_intra_filter['Protein2']
            PPI_library_intra_filter = PPI_library_intra_filter.drop(['Site1', 'Site2'], axis=1)
            pair_counts_intra = PPI_library_intra_filter['Pair'].value_counts().reset_index()
            pair_counts_intra.columns = ['Pair', 'Frequency']
        
        # 处理 inter 数据
        if not data_inter.empty:
            PPI_library_inter_filter = PPI_library_inter.str.extract(pattern)
            PPI_library_inter_filter.columns = ['Protein1', 'Site1', 'Protein2', 'Site2']
            mask = PPI_library_inter_filter['Protein1'] > PPI_library_inter_filter['Protein2']
            PPI_library_inter_filter.loc[mask, ['Protein1', 'Protein2']] = PPI_library_inter_filter.loc[mask, ['Protein2', 'Protein1']].values
            PPI_library_inter_filter.loc[mask, ['Site1', 'Site2']] = PPI_library_inter_filter.loc[mask, ['Site2', 'Site1']].values
            PPI_library_inter_filter['Pair'] = PPI_library_inter_filter['Protein1'] + '-' + PPI_library_inter_filter['Protein2']
            PPI_library_inter_filter = PPI_library_inter_filter.drop(['Site1', 'Site2'], axis=1)
            pair_counts_inter = PPI_library_inter_filter['Pair'].value_counts().reset_index()
            pair_counts_inter.columns = ['Pair', 'Frequency']

       # 处理 intra peptide 数据
        if not data_intra.empty:
            data_peptide_intra = data_intra[['Peptide', 'Modifications', 'Proteins']].copy().reset_index(drop=True)
            all_peptide_info = []
        
            for i in range(len(data_peptide_intra)):
                dfi = data_peptide_intra['Proteins'].iloc[i].split('/')
                dfi = pd.Series(dfi).str.strip()
                dfi = dfi[dfi != '']
                dfi = dfi.str.extract(pattern)
                dfi.columns = ['Protein1', 'Site1', 'Protein2', 'Site2']
                dfi = dfi[dfi['Protein1'] == dfi['Protein2']]
                dfi['Pair'] = dfi['Protein1'] + '-' + dfi['Protein2']
                dfi = pd.merge(dfi, pair_counts_intra, on='Pair', how='left')
                dfi['Frequency'] = dfi['Frequency'].fillna(0)
                if not dfi.empty:
                    for _, row in dfi.iterrows():
                        all_peptide_info.append([data_peptide_intra['Peptide'].iloc[i],
                                                 data_peptide_intra['Modifications'].iloc[i],
                                                 data_peptide_intra['Proteins'].iloc[i],
                                                 row['Protein1'], row['Protein2'],
                                                 int(row['Site1']), int(row['Site2']),
                                                 row['Frequency']])
        
            data_peptide_intra = pd.DataFrame(all_peptide_info, columns=['Peptide', 'Modifications', 'Proteins', 'Protein1', 'Protein2', 'Site1', 'Site2', 'Frequency'])

        

        # 处理 inter peptide 数据
        if not data_inter.empty:
            data_peptide_inter = data_inter[['Peptide', 'Modifications', 'Proteins']].copy().reset_index(drop=True)
            data_peptide_inter[['Protein1', 'Protein2', 'Site1', 'Site2', 'Frequency']] = ''
            for j in range(len(data_peptide_inter)):
                dfj = data_peptide_inter['Proteins'].iloc[j].split('/')
                dfj = pd.Series(dfj).str.strip()
                dfj = dfj[dfj != '']
                dfj = dfj.str.extract(pattern)
                dfj.columns = ['Protein1', 'Site1', 'Protein2', 'Site2']
                mask1 = dfj['Protein1'] > dfj['Protein2']
                dfj.loc[mask1, ['Protein1', 'Protein2']] = dfj.loc[mask1, ['Protein2', 'Protein1']].values
                dfj.loc[mask1, ['Site1', 'Site2']] = dfj.loc[mask1, ['Site2', 'Site1']].values
                dfj['Pair'] = dfj['Protein1'] + '-' + dfj['Protein2']
                dfj = pd.merge(dfj, pair_counts_inter, on='Pair', how='left')
                dfj['Frequency'] = dfj['Frequency'].fillna(0)
                dfj = dfj.sort_values(by='Frequency', ascending=False)
                if not dfj.empty:
                    data_peptide_inter.loc[j, 'Protein1'] = dfj.iloc[0]['Protein1']
                    data_peptide_inter.loc[j, 'Protein2'] = dfj.iloc[0]['Protein2']
                    data_peptide_inter.loc[j, 'Site1'] = int(dfj.iloc[0]['Site1'])
                    data_peptide_inter.loc[j, 'Site2'] = int(dfj.iloc[0]['Site2'])
                    data_peptide_inter.loc[j, 'Frequency'] = dfj.iloc[0]['Frequency']

        # 过滤并处理 data_peptide_inter 和 data_site_inter
        if not data_inter.empty:
            data_peptide_inter = data_peptide_inter[data_peptide_inter['Protein1'] != '']
            data_site_inter = data_peptide_inter.groupby(['Protein1', 'Protein2', 'Site1', 'Site2']).size().reset_index(name='Spectra_Count')
            all_data_peptide_inter.append(data_peptide_inter)
            all_data_site_inter.append(data_site_inter)
        
        # 过滤并处理 data_peptide_intra 和 data_site_intra
        if not data_intra.empty:
            data_site_intra = data_peptide_intra.groupby(['Protein1', 'Protein2', 'Site1', 'Site2']).size().reset_index(name='Spectra_Count')
            data_site_intra = data_site_intra[data_site_intra['Protein1'] != '']
            all_data_peptide_intra.append(data_peptide_intra)
            all_data_site_intra.append(data_site_intra)       
# 合并所有文件的data_peptide、 data_site 和 data_ppi
merged_data_peptide_intra = pd.concat(all_data_peptide_intra, ignore_index=True)
merged_data_peptide_inter = pd.concat(all_data_peptide_inter, ignore_index=True)
merged_data_site_intra = pd.concat(all_data_site_intra, ignore_index=True)
merged_data_site_inter = pd.concat(all_data_site_inter, ignore_index=True)
# 在合并后的 data_site 上进行去重和 'Spectra_Count' 求和
merged_data_site_intra = merged_data_site_intra.groupby(['Protein1', 'Protein2', 'Site1', 'Site2'])['Spectra_Count'].sum().reset_index()
merged_data_ppi_intra = merged_data_site_intra.copy()
merged_data_ppi_intra['Site_Count'] = merged_data_site_intra.groupby(['Protein1', 'Protein2'])['Protein1'].transform('count')
merged_data_ppi_intra = merged_data_ppi_intra.groupby(['Protein1', 'Protein2']).agg({'Spectra_Count': 'sum', 'Site_Count': 'first'}).reset_index()

merged_data_site_inter = merged_data_site_inter.groupby(['Protein1', 'Protein2', 'Site1', 'Site2'])['Spectra_Count'].sum().reset_index()
merged_data_ppi_inter = merged_data_site_inter.copy()
merged_data_ppi_inter['Site_Count'] = merged_data_site_inter.groupby(['Protein1', 'Protein2'])['Protein1'].transform('count')
merged_data_ppi_inter = merged_data_ppi_inter.groupby(['Protein1', 'Protein2']).agg({'Spectra_Count': 'sum', 'Site_Count': 'first'}).reset_index()

# 复制并去重
unique_peptide_intra = merged_data_peptide_intra[['Peptide', 'Modifications']].drop_duplicates()
unique_peptide_inter = merged_data_peptide_inter[['Peptide', 'Modifications']].drop_duplicates()

# 假设 uniprot-gene 文件路径
uniprot_gene_file = "/mnt/Stringent filtering of cross-links/data/uniprot-gene-20220927.xlsx"

# 加载 uniprot-gene 数据
uniprot_gene_df = pd.read_excel(uniprot_gene_file)
# 提取第一个基因名称
uniprot_gene_df['Gene_name'] = uniprot_gene_df['Gene'].str.split().str[0]
# 创建映射字典
uniprot_to_gene = uniprot_gene_df.set_index('Entry')['Gene_name'].to_dict()

# 添加 Gene1 和 Gene2 列的函数
def add_gene_columns(df):
    df['Gene1'] = df['Protein1'].map(uniprot_to_gene)
    df['Gene2'] = df['Protein2'].map(uniprot_to_gene)
    return df

# 为 intra 和 inter 数据表添加 Gene 列
merged_data_site_intra = add_gene_columns(merged_data_site_intra)
merged_data_ppi_intra = add_gene_columns(merged_data_ppi_intra)
merged_data_site_inter = add_gene_columns(merged_data_site_inter)
merged_data_ppi_inter = add_gene_columns(merged_data_ppi_inter)

explanations = {
    'intra_protein_peptides': (
        "Identified non-redundant intra protein crosslink peptides"
    ),
    'intra_protein_peptides_to_sites': (
        "Results of matching intra protein crosslink peptides to intra protein site pairs, "
        "Frequency indicates the identification frequency (number of spectra) of protein in all identifications"
    ),
    'intra_protein_sites': (
        "Identified non-redundant intra protein site pairs, "
        "Spectra_Count indicates the number of spectra corresponding to the site pairs"
    ),
    'intra_proteins': (
        "Identified non-redundant intra proteins, "
        "Spectra_Count indicates the number of spectra corresponding to intra proteins, "
        "Site_Count indicates the number of intra protein site pairs identified for each intra protein"
    ),
    'inter_protein_peptides': (
        "Identified non-redundant inter protein crosslink peptides"
    ),
    'inter_protein_peptides_to_sites': (
        "Results of matching inter protein crosslink peptides to inter protein site pairs, "
        "Frequency indicates the identification frequency (number of spectra) of PPI in all identifications. "
        "Retain only the inter protein site pair with the highest identification frequency for each pair of cross-linked peptide"
    ),
    'inter_protein_sites': (
        "Identified non-redundant inter protein site pairs, "
        "Spectra_Count indicates the number of spectra corresponding to the site pairs"
    ),
    'inter_protein_PPI': (
        "Identified non-redundant PPI, "
        "Spectra_Count indicates the number of spectra corresponding to PPI, "
        "Site_Count indicates the number of inter protein site pairs identified for each PPI"
    )
}

# 将所有数据写入一个Excel文件的不同Sheet中
with pd.ExcelWriter('Tau_plink(inter_high_confidence).xlsx') as writer:
    for sheet_name, explanation in explanations.items():
        if sheet_name == 'intra_protein_peptides':
            df = unique_peptide_intra
        elif sheet_name == 'intra_protein_peptides_to_sites':
            df = merged_data_peptide_intra
        elif sheet_name == 'intra_protein_sites':
            df = merged_data_site_intra
        elif sheet_name == 'intra_proteins':
            df = merged_data_ppi_intra
        elif sheet_name == 'inter_protein_peptides':
            df = unique_peptide_inter
        elif sheet_name == 'inter_protein_peptides_to_sites':
            df = merged_data_peptide_inter
        elif sheet_name == 'inter_protein_sites':
            df = merged_data_site_inter
        elif sheet_name == 'inter_protein_PPI':
            df = merged_data_ppi_inter

        # 创建包含解释说明的DataFrame
        explanation_df = pd.DataFrame([explanation], columns=[df.columns[0]])
        # 将解释说明和数据写入Excel文件
        explanation_df.to_excel(writer, sheet_name=sheet_name, index=False, header=False)
        df.to_excel(writer, sheet_name=sheet_name, index=False, startrow=1)

