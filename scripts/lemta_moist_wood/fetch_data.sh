#!/usr/bin/env bash
# Download the LEMTA moist-wood dataset (doi:10.57745/CBLZ8P) file by file
# through the Recherche Data Gouv API into $LEMTA_DATA (default
# ~/Documents/data/lemta-CBLZ8P) and check each file against the MD5 the
# repository states, as pinned in data_files.toml. Files already present with
# the right checksum are kept.
#
#   bash scripts/lemta_moist_wood/fetch_data.sh
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
manifest="$here/data_files.toml"
dest="${LEMTA_DATA:-$HOME/Documents/data/lemta-CBLZ8P}"
api="$(sed -n 's/^api = "\(.*\)"$/\1/p' "$manifest")"
mkdir -p "$dest"

# name, id and md5 of every [[files]] entry, one tab-separated line each
entries="$(awk -F' = ' '
    /^\[\[files\]\]/ { if (name != "") print name "\t" id "\t" md5; name = id = md5 = "" }
    $1 == "name" { name = $2; gsub(/"/, "", name) }
    $1 == "id"   { id = $2 }
    $1 == "md5"  { md5 = $2; gsub(/"/, "", md5) }
    END { if (name != "") print name "\t" id "\t" md5 }
' "$manifest")"

failed=0
while IFS=$'\t' read -r name id md5; do
    path="$dest/$name"
    if [[ -f "$path" ]] && [[ "$(md5sum "$path" | cut -d' ' -f1)" == "$md5" ]]; then
        continue
    fi
    echo "fetching $name"
    curl -sSL --fail -m 300 -o "$path.part" "$api/access/datafile/$id"
    got="$(md5sum "$path.part" | cut -d' ' -f1)"
    if [[ "$got" == "$md5" ]]; then
        mv "$path.part" "$path"
    else
        echo "checksum mismatch for $name: got $got, expected $md5" >&2
        rm -f "$path.part"
        failed=1
    fi
done <<< "$entries"

if [[ $failed -ne 0 ]]; then
    echo "some files failed their checksum" >&2
    exit 1
fi
echo "dataset complete in $dest"
