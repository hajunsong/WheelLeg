#!/usr/bin/env bash
# Linux / WSL / macOS 공용 셋업 스크립트.
# Windows(PowerShell) 는 README.md 참고.
#
#   ./setup.sh              # python3 사용
#   PYTHON=python3.12 ./setup.sh
set -euo pipefail
cd "$(dirname "$0")"

PY=${PYTHON:-python3}
VENV=.venv-linux

command -v "$PY" >/dev/null 2>&1 || {
    echo "오류: '$PY' 를 찾을 수 없습니다." >&2; exit 1; }

"$PY" -c 'import venv' >/dev/null 2>&1 || {
    echo "오류: venv 모듈이 없습니다." >&2
    echo "  Debian/Ubuntu : sudo apt install python3-venv" >&2
    echo "  Fedora/RHEL   : sudo dnf install python3-virtualenv" >&2
    exit 1; }

echo "[1/2] $VENV 생성 ($("$PY" --version))"
rm -rf "$VENV"
"$PY" -m venv "$VENV"

echo "[2/2] 의존성 설치"
"$VENV/bin/python" -m pip install --quiet --upgrade pip
"$VENV/bin/python" -m pip install --quiet -r requirements.txt

echo
"$VENV/bin/python" -c 'import sys, numpy, matplotlib; print(f"  Python {sys.version.split()[0]} / numpy {numpy.__version__} / matplotlib {matplotlib.__version__}")'
echo
echo "완료.  실행 :  ./$VENV/bin/python main.py"
echo "  DISPLAY 가 없으면 compare.png / error.png 로 저장됩니다."
