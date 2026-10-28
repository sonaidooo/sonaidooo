#!/bin/bash
set -e

samples="data/reference_data/samples.txt"
out_dir="data/1000_genomes_data"
temp_dir="${out_dir}/temp_chunks"

mkdir -p "$out_dir" "$temp_dir"

chr_lens=(
  0
  249000000 243000000 199000000 191000000 182000000
  172000000 160000000 146000000 139000000 134000000
  136000000 134000000 115000000 108000000 102000000
  91000000  84000000  81000000  59000000  65000000
  47000000  51000000
)

render_overall_bar() {
  local current=$1
  local total=$2
  local chr=$3
  local width=20
  local percent=$(( current * 100 / total ))
  local filled=$(( current * width / total ))
  local empty=$(( width - filled ))

  local bar=""
  for ((j=0; j<filled; j++)); do bar="${bar}█"; done
  for ((j=0; j<empty; j++)); do bar="${bar}░"; done

  printf "Chr %-2s Overall: [%s] %2d/%2d chunks (%3d%%)\n" "$chr" "$bar" "$current" "$total" "$percent"
}

for i in {2..22}; do
  #resume in case of dropped connection
  if [ -f "${out_dir}/subset_chr${i}.vcf.gz.tbi" ]; then
    echo "Chromosome ${i} already complete. Skipping..."
    continue
  fi

  echo -e "\n=========================================="
  echo "Processing Chromosome ${i}..."
  echo "=========================================="
  
  url="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.recalibrated_variants.vcf.gz"
  chr_max=${chr_lens[$i]}
  
  total_chunks=0
  for ((start=1; start<=chr_max; start+=20000000)); do
    total_chunks=$((total_chunks + 1))
  done

  chunk_list=()
  chunk_index=1

  for ((start=1; start<=chr_max; start+=20000000)); do
    end=$((start + 19999999))
    if [ "$end" -gt "$chr_max" ]; then end=$chr_max; fi
    
    region="chr${i}:${start}-${end}"
    chunk_file="${temp_dir}/chr${i}_chunk_${chunk_index}.vcf.gz"

    if [ -f "${chunk_file}.done" ]; then
      render_overall_bar "$chunk_index" "$total_chunks" "$i"
      chunk_list+=("$chunk_file")
      chunk_index=$((chunk_index + 1))
      continue
    fi

    render_overall_bar "$((chunk_index-1))" "$total_chunks" "$i"

    success=0
    for attempt in {1..5}; do
      if bcftools view -r "$region" -S "$samples" --force-samples -Ou "$url" 2>/dev/null | \
         bcftools annotate -x FORMAT/AD,FORMAT/DP,FORMAT/GQ,FORMAT/PL,FORMAT/MIN_DP,FORMAT/SB -Oz | \
         python3 -c '
import sys, time
c_idx, t_chunks, reg = sys.argv[1], sys.argv[2], sys.argv[3]
t0 = time.time()
nbytes = 0
est = 20 * 1024 * 1024
while True:
    b = sys.stdin.buffer.read(65536)
    if not b: break
    sys.stdout.buffer.write(b)
    nbytes += len(b)
    dt = time.time() - t0
    mb = nbytes / 1048576
    speed = mb / dt if dt > 0 else 0
    pct = min(int((nbytes / est) * 100), 99)
    filled = pct // 5
    bar = "█" * filled + "░" * (20 - filled)
    sys.stderr.write("\r  └─ Chunk %s/%s [%s]: [%s] %.1f MB (%.2f MB/s)" % (c_idx, t_chunks, reg, bar, mb, speed))
    sys.stderr.flush()
mb = nbytes / 1048576
sys.stderr.write("\r  └─ Chunk %s/%s [%s]: [████████████████████] %.1f MB (Done!)              \n" % (c_idx, t_chunks, reg, mb))
' "$CHUNK_IDX" "$TOTAL_CHUNKS" "$REGION" > "$CHUNK_FILE"; then
        
        if [ -s "$chunk_file" ]; then
          bcftools index -t "$chunk_file" 2>/dev/null || true
          touch "${chunk_file}.done"
          success=1
          break
        fi
      fi
      
      echo -e "\n  [!] Network drop on ${region} (Attempt ${attempt}/5). Retrying in 3s..."
      rm -f "$chunk_file"*
      sleep 3
    done

    if [ $success -eq 0 ]; then
      echo -e "\nError: Failed chunk ${region} after 5 attempts."
      exit 1
    fi

    chunk_list+=("$chunk_file")
    chunk_index=$((chunk_index + 1))
  done

  render_overall_bar "$total_chunks" "$total_chunlks" "$i"
  echo "Concatenating chunks for Chromosome ${i}..."
  bcftools concat -a -Oz -o "${out_dir}/subset_chr${i}.vcf.gz" "${chunk_list[@]}"
  
  echo "Indexing subset_chr${i}.vcf.gz..."
  bcftools index -t "${out_dir}/subset_chr${i}.vcf.gz"

  rm -f "${temp_dir}/chr${i}_chunk_"*
  echo "Chromosome ${i} complete!"
done

rm -rf "$temp_dir"
echo "All 22 autosomes finished successfully!"