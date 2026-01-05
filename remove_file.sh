#!/usr/bin/env bash

# 简化修复版：增强输入检查、处理 CRLF、在删除失败时报错
set -u

if [ $# -eq 0 ]; then
    echo "Usage: $0 <filename_list>" >&2
    exit 1
fi

# 检查输入文件是否存在且可读
if [ ! -r "$1" ]; then
    echo "Error: List file '$1' not found or not readable." >&2
    exit 1
fi

# 读取文件中的每一行并处理；保留前导空格但跳过空白行；去掉可能的 CR (\r)
while IFS= read -r file || [ -n "$file" ]; do
    # 删除行尾可能存在的 CR（Windows CRLF 行尾）
    file="${file//$'\r'/}"

    # 跳过空行或仅包含空白字符的行
    if [ -z "$(printf "%s" "$file" | tr -d '[:space:]')" ]; then
        continue
    fi

    # 检查路径是否为空（再次保险）
    if [ -z "$file" ]; then
        continue
    fi

    # 检查文件或目录是否存在
    if [ -e "$file" ]; then
        rm -- "$file" 2>/dev/null
        if [ $? -eq 0 ]; then
            echo "Removed: $file"
        else
            echo "Failed to remove: $file" >&2
        fi
    else
        echo "Not found, skipping: $file"
    fi
done < "$1"

echo "File removal process completed."