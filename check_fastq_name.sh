#!/bin/bash

# 方法3：简单检查
echo "查找二级编码相同但一级编码不同且三级编码相同的文件..."
echo "================================================"

# 使用关联数组跟踪：键为"二级编码_三级编码"，值为一级编码
declare -A seen_patterns

for file in *_*_*_*.fq.gz; do
    if [[ -f "$file" ]]; then
        # 提取一级编码（前两个字段：E251211007_L01）
        primary=$(echo "$file" | cut -d'_' -f1-2)
        
        # 提取二级编码（第三个字段：DY-99）
        secondary=$(echo "$file" | cut -d'_' -f3)
        
        # 提取三级编码（第四个字段去掉.fq.gz：2）
        tertiary=$(echo "$file" | cut -d'_' -f4 | sed 's/\.fq\.gz$//')
        
        # 构造复合键：二级编码_三级编码
        composite_key="${secondary}_${tertiary}"
        
        # 检查这个复合键是否已经出现过
        if [[ -n "${seen_patterns[$composite_key]}" ]]; then
            # 如果已经出现过，获取之前记录的一级编码
            existing_primary="${seen_patterns[$composite_key]}"
            
            # 如果一级编码不同，则视为冲突
            if [[ "$existing_primary" != "$primary" ]]; then
                echo "冲突: 二级编码 $secondary 和三级编码 $tertiary"
                echo "  文件1: ${existing_primary}_${secondary}_${tertiary}.fq.gz（模式）"
                echo "  文件2: $file"
                echo ""
            fi
            # 如果一级编码相同，不视为冲突
        else
            # 第一次见到这个复合键，记录它的一级编码
            seen_patterns["$composite_key"]="$primary"
        fi
    fi
done