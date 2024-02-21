'''
Author: ShengyiXu xushengyichn@outlook.com
Date: 2024-02-21 14:21:23
LastEditors: ShengyiXu xushengyichn@outlook.com
LastEditTime: 2024-02-21 14:56:17
FilePath: \manuscriptf:\git\xihoumen_inverse_force_estimation\20231203 third edition\remove_pic.py
Description: 

Copyright (c) 2024 by ${git_name_email}, All Rights Reserved. 
'''

# %%
import os
import re
import shutil
from bs4 import BeautifulSoup

# 检查delete文件夹中是否有文件的函数
def check_and_delete_files_in_delete_folder(delete_folder_path):
    if os.path.exists(delete_folder_path) and os.listdir(delete_folder_path):
        # 如果delete文件夹存在且不为空，则提示用户
        response = input("delete文件夹中有文件，是否要删除这些文件？(y/n): ").strip().lower()
        if response == 'y':
            for file in os.listdir(delete_folder_path):
                file_path = os.path.join(delete_folder_path, file)
                try:
                    if os.path.isfile(file_path):  # 确保是文件
                        os.remove(file_path)
                    elif os.path.isdir(file_path):  # 如果有子目录，可以选择递归删除
                        shutil.rmtree(file_path)
                    print(f"已删除：{file}")
                except Exception as e:
                    print(f"删除文件{file}失败：{e}")
            print("已清空delete文件夹中的所有文件。")
        else:
            print("保留了delete文件夹中的所有文件。")
    else:
        print("delete文件夹为空或不存在，无需进行删除操作。")
        
def process_markdown_file(md_file_path, working_directory):
    md_name_without_ext = os.path.splitext(os.path.basename(md_file_path))[0]
    assets_directory_name = md_name_without_ext + '.assets'
    assets_directory_path = os.path.join(working_directory, assets_directory_name)
    delete_path = os.path.join(assets_directory_path, 'delete')

    if not os.path.exists(assets_directory_path):
        os.makedirs(assets_directory_path)
        print(f"已创建.assets目录：{assets_directory_path}")

    if not os.path.exists(delete_path):
        os.makedirs(delete_path)
        print(f"已创建delete目录：{delete_path}")

    pattern = re.compile(r'!\[.*]\((.*)\)')
    references = set()

    with open(md_file_path, 'r', encoding='utf-8') as f:
        text = f.read()
        references.update(pattern.findall(text))
        soup = BeautifulSoup(text, 'html.parser')
        for img in soup.find_all('img'):
            src = img.get('src')
            if src:
                references.add(os.path.join(os.path.dirname(md_file_path), src))

    usedImages = {os.path.basename(each) for each in references}

    allImageNum = 0
    deleteNum = 0
    for file in os.listdir(assets_directory_path):
        file_path = os.path.join(assets_directory_path, file)
        if os.path.isfile(file_path):  # 确保是文件而不是目录
            allImageNum += 1
            if file not in usedImages and file != 'delete':  # 排除delete目录
                shutil.move(file_path, os.path.join(delete_path, file))
                print(f"移动图片{file}到delete文件夹中")
                deleteNum += 1


    print(f"{md_file_path}中：.assets文件夹中总共有图片{allImageNum}张，总共移动{deleteNum}张图片到delete中")
    delete_folder_path = os.path.join(assets_directory_path, 'delete')
    check_and_delete_files_in_delete_folder(delete_folder_path)
    

working_directory = input("请输入Markdown文件所在的目录路径（留空使用当前目录）: ").strip() or os.getcwd()
md_files = [file for file in os.listdir(working_directory) if file.endswith('.md')]

if not md_files:
    print("没有找到Markdown文件。")
else:
    print("找到以下Markdown文件：")
    for i, file in enumerate(md_files, start=1):
        print(f"{i}. {file}")
    print(f"{len(md_files) + 1}. 处理所有Markdown文件")

    choice = int(input("请输入您想操作的选项编号："))

    if choice == len(md_files) + 1:
        for md_file in md_files:
            md_file_path = os.path.join(working_directory, md_file)
            process_markdown_file(md_file_path, working_directory)
    elif 0 < choice <= len(md_files):
        selected_md_file = md_files[choice - 1]
        md_file_path = os.path.join(working_directory, selected_md_file)
        process_markdown_file(md_file_path, working_directory)
    else:
        print("无效的选择。")

# %%
