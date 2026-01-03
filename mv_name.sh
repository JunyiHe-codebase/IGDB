#!/bin/bash

# 遍历当前目录下所有.fq.gz文件
for file in *_*_*_*.fq.gz; do
    # 检查文件是否存在（避免空匹配）
    if [[ -f "$file" ]]; then
        # 提取二级编码（倒数第二个下划线后的部分）
        # 例如：从E251028002_L01_DY-109_1.fq.gz中提取DY-109
        secondary_code=$(echo "$file" | grep -o '_DY-[^_]*' | cut -d'_' -f2)
        
        # 提取三级编码（最后一个下划线后的部分，去除.fq.gz）
        # 例如：从E251028002_L01_DY-109_1.fq.gz中提取1
        tertiary_code=$(echo "$file" | sed 's/.*_\([^.]*\)\.fq\.gz/\1/')
        
        # 构建新文件名
        new_name="${secondary_code}_${tertiary_code}.fq.gz"
        
        # 检查新文件名是否已经存在
        if [[ -e "$new_name" ]]; then
            echo "跳过: $file -> $new_name (目标文件已存在)"
        else
            # 重命名文件
            mv "$file" "$new_name"
            echo "已重命名: $file -> $new_name"
        fi
    fi
done

echo "批量重命名完成！"