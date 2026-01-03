#!/bin/bash

for file in *.zip; do
    # 检查是否存在zip文件
    [ -f "$file" ] || continue
    
    # 解压zip文件
    echo "正在处理: $file"
    unzip -q "$file"

    dir_name="${file%.zip}"
    
    # 检查目录是否存在
    if [ -d "$dir_name" ]; then
        report_file="$dir_name/summary.txt"
        
        if [ -f "$report_file" ]; then
            # 在报告中查找相关信息，并标记来源
            echo "=== 来自: $file ===" >> total_summary.txt
            grep -A 5 "Overrepresented sequences" "$report_file" >> total_summary.txt
            grep -A 5 "Adapter Content" "$report_file" >> total_summary.txt
            echo "" >> total_summary.txt  # 添加空行分隔
        else
            echo "警告: 在 $dir_name 中未找到 summary.txt"
        fi
    else
        echo "警告: 解压后未找到预期的目录: $dir_name"
    fi
done