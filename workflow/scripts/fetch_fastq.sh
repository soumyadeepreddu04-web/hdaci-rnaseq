#!/usr/bin/env bash
# Fetch one paired-end library from ENA.
#
#   fetch_fastq.sh URL_1 URL_2 OUT_1 OUT_2 PAIRS WHOLE_FILE MD5_1 MD5_2
#
# WHOLE_FILE = True   download both files completely and check them against the ENA md5 sums.
# WHOLE_FILE = False  stream each file and keep the first PAIRS records. Mates are stored in the
#                     same order in the two files, so the first N records of each stay paired.
set -u

url1=$1; url2=$2; out1=$3; out2=$4; pairs=$5; whole_file=$6; md5_1=$7; md5_2=$8
lines=$((pairs * 4))

checksum() {
    if command -v md5sum > /dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi
}

if [ "$whole_file" = "True" ]; then
    curl -sS -L --retry 10 --retry-delay 15 -C - -o "$out1" "$url1" || exit 1
    curl -sS -L --retry 10 --retry-delay 15 -C - -o "$out2" "$url2" || exit 1
    [ "$(checksum "$out1")" = "$md5_1" ] || { echo "md5 mismatch: $out1"; exit 1; }
    [ "$(checksum "$out2")" = "$md5_2" ] || { echo "md5 mismatch: $out2"; exit 1; }
    exit 0
fi

# head closes the pipe once it has enough lines, which makes curl exit non-zero by design,
# so success is judged from the number of lines written rather than from the exit status.
# ENA connections stall or drop now and then, hence the stall timeout and the retries.
stream() {
    local url=$1 out=$2 attempt n
    for attempt in 1 2 3 4 5 6; do
        curl -sS -L --connect-timeout 30 --speed-limit 10000 --speed-time 60 "$url" 2> /dev/null \
            | gzip -dc 2> /dev/null | head -n "$lines" | gzip -1 > "$out"
        n=$(gzip -dc "$out" | wc -l | tr -d ' ')
        echo "$(basename "$out"): attempt $attempt, $n of $lines lines"
        [ "$n" -eq "$lines" ] && return 0
        sleep 60
    done
    return 1
}

stream "$url1" "$out1" & pid1=$!
stream "$url2" "$out2" & pid2=$!
status=0
wait "$pid1" || status=1
wait "$pid2" || status=1
[ "$status" -eq 0 ] || echo "incomplete download"
exit "$status"
