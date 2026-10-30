#!/usr/bin/env bash
# chunk_fetch.sh URL BYTES OUTFILE
# Download one file in fixed-size range chunks (default 128 MB; CHUNK_MB, or CHUNK_BYTES for tests).
# A chunk is accepted ONLY if the server answers HTTP 206 and exactly the expected number of bytes arrive;
# otherwise the chunk is discarded and fetched again from its own start (no appending to partial data, so no splices).
# Finished chunks stay in OUTFILE.chunks/ : an interrupted run continues at chunk level. At the end the chunks are
# joined in order into OUTFILE and the size is checked (the caller checks the md5).
# After HTTP 403/429 the wait is 180 s, otherwise 15 s. Up to CHUNK_TRIES (default 15) tries per chunk.
set -uo pipefail
URL=${1:?url}; BYTES=${2:?bytes}; OUTF=${3:?outfile}
CHUNK=${CHUNK_BYTES:-$(( ${CHUNK_MB:-128} * 1048576 ))}; TRIES=${CHUNK_TRIES:-15}
D="$OUTF.chunks"; mkdir -p "$D"
N=$(( (BYTES + CHUNK - 1) / CHUNK ))
for ((i=0; i<N; i++)); do
  s=$(( i*CHUNK )); e=$(( s+CHUNK-1 )); [ "$e" -ge "$BYTES" ] && e=$(( BYTES-1 ))
  want=$(( e-s+1 )); cf=$(printf '%s/%05d' "$D" "$i")
  [ -f "$cf" ] && [ "$(stat -c %s "$cf")" = "$want" ] && continue
  ok=0
  for ((t=1; t<=TRIES; t++)); do
    rm -f "$cf.tmp"
    code=$(curl -sSL --connect-timeout 30 --speed-time 60 --speed-limit 20000 -r "$s-$e" -o "$cf.tmp" -w '%{http_code}' "$URL" 2>/dev/null) || code=000
    got=$(stat -c %s "$cf.tmp" 2>/dev/null || echo 0)
    if [ "$code" = 206 ] && [ "$got" = "$want" ]; then mv "$cf.tmp" "$cf"; ok=1; break; fi
    rm -f "$cf.tmp"
    w=${WAIT_SHORT:-15}; { [ "$code" = 403 ] || [ "$code" = 429 ]; } && w=${WAIT_LONG:-180}
    echo "$(date -u +%T) chunk $((i+1))/$N try $t failed (http $code, got $got of $want bytes); waiting $w s" >&2
    sleep "$w"
  done
  [ "$ok" = 1 ] || { echo "chunk $((i+1))/$N failed after $TRIES tries (chunks kept; re-run continues)" >&2; exit 1; }
  { [ $(( (i+1) % 8 )) = 0 ] || [ $((i+1)) = "$N" ]; } && echo "$(date -u +%T) chunk $((i+1))/$N done" >&2
done
rm -f "$OUTF.assembling"
cat "$D"/[0-9][0-9][0-9][0-9][0-9] > "$OUTF.assembling" || exit 1
[ "$(stat -c %s "$OUTF.assembling")" = "$BYTES" ] || { echo "assembled size wrong" >&2; rm -f "$OUTF.assembling"; exit 1; }
mv -f "$OUTF.assembling" "$OUTF" && rm -rf "$D"
