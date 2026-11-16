#!/usr/bin/env bash
# Download the LEMTA dry argon dataset (doi:10.57745/QRYAJ9) file by file
# through the Recherche Data Gouv API into $LEMTA_ARGON_DATA (default
# ~/Documents/data/lemta-QRYAJ9). Each original upload is checked against the
# MD5 the repository states; each spreadsheet's tab-delimited export (the
# repository's default download, which the loader reads) is checked against
# the SHA-256 pinned in data_files.toml. Files already present with the right
# checksum are kept.
#
#   bash scripts/lemta_argon/fetch_data.sh
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
manifest="$here/data_files.toml"
dest="${LEMTA_ARGON_DATA:-$HOME/Documents/data/lemta-QRYAJ9}"
api="$(sed -n 's/^api = "\(.*\)"$/\1/p' "$manifest")"
mkdir -p "$dest"

# name, id, md5, tab and tab_sha256 of every [[files]] entry, one
# tab-separated line each ("-" where an entry has no tab export)
entries="$(awk -F' = ' '
    function flush() { if (name != "") print name "\t" id "\t" md5 "\t" tab "\t" tsha }
    /^\[\[files\]\]/ { flush(); name = id = md5 = ""; tab = tsha = "-" }
    $1 == "name"       { name = $2; gsub(/"/, "", name) }
    $1 == "id"         { id = $2 }
    $1 == "md5"        { md5 = $2; gsub(/"/, "", md5) }
    $1 == "tab"        { tab = $2; gsub(/"/, "", tab) }
    $1 == "tab_sha256" { tsha = $2; gsub(/"/, "", tsha) }
    END { flush() }
' "$manifest")"

failed=0
# fetch <url> <path> <hash command> <expected>
fetch() {
    local url="$1" path="$2" hash="$3" want="$4"
    if [[ -f "$path" ]] && [[ "$($hash "$path" | cut -d' ' -f1)" == "$want" ]]; then
        return 0
    fi
    echo "fetching $(basename "$path")"
    curl -sSL --fail -m 300 -o "$path.part" "$url"
    local got
    got="$($hash "$path.part" | cut -d' ' -f1)"
    if [[ "$got" == "$want" ]]; then
        mv "$path.part" "$path"
    else
        echo "checksum mismatch for $(basename "$path"): got $got, expected $want" >&2
        rm -f "$path.part"
        failed=1
    fi
}

while IFS=$'\t' read -r name id md5 tab tsha; do
    if [[ "$tab" == "-" ]]; then
        fetch "$api/access/datafile/$id" "$dest/$name" md5sum "$md5"
    else
        fetch "$api/access/datafile/$id?format=original" "$dest/$name" md5sum "$md5"
        fetch "$api/access/datafile/$id" "$dest/$tab" sha256sum "$tsha"
    fi
done <<< "$entries"

if [[ $failed -ne 0 ]]; then
    echo "some files failed their checksum" >&2
    exit 1
fi
echo "dataset complete in $dest"
