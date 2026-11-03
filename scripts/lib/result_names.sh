#!/usr/bin/env bash
# Shared validation for user-visible run directories, not internal task folders.
set -euo pipefail

require_date_first_result_name() {
    local result_name="${1:-}" legacy_mode="${2:-}" result_date normalized_date
    if [[ "$result_name" =~ ^([0-9]{4}-[0-9]{2}-[0-9]{2})_([01][0-9]|2[0-3])-([0-5][0-9])-[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        result_date="${BASH_REMATCH[1]}"
    elif [[ "$legacy_mode" == legacy-existing && "$result_name" =~ ^([0-9]{4}-[0-9]{2}-[0-9]{2})[-_][A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        result_date="${BASH_REMATCH[1]}"
    elif [[ "$legacy_mode" == legacy-existing && "$result_name" =~ ^([0-9]{4})([0-9]{2})([0-9]{2})[-_][A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        result_date="${BASH_REMATCH[1]}-${BASH_REMATCH[2]}-${BASH_REMATCH[3]}"
    else
        echo "New result name must start with UTC date and time: YYYY-MM-DD_HH-MM-description-vN: $result_name" >&2
        return 1
    fi
    if ! normalized_date="$(date -d "$result_date" +%F 2>/dev/null)" || [[ "$normalized_date" != "$result_date" ]]; then
        echo "Invalid calendar date in result name: $result_name" >&2
        return 1
    fi
}
