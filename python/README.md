# python — matlab/ 의 1:1 이식

MATLAB 코드와 파일명·함수명·변수명을 그대로 유지해 옮긴 것입니다.

## 실행

Windows 용(`.venv`)과 WSL/리눅스 용(`.venv-wsl`) venv 가 각각 들어 있습니다.
`.venv` 는 `Scripts/*.exe`, `.venv-wsl` 은 `bin/` 구조라 서로 호환되지 않아 분리했습니다.

### Windows (PowerShell)

```powershell
cd C:\Users\gkwns\Desktop\WheelLeg\python
.\.venv\Scripts\python.exe main.py
```

활성화해서 쓰려면 `.\.venv\Scripts\Activate.ps1` 실행 후 `python main.py`, 나올 때 `deactivate`.
VSCode 에서는 `Ctrl+Shift+P` → **Python: Select Interpreter** 로
`.\python\.venv\Scripts\python.exe` 를 지정하면 실행 버튼이 그대로 동작합니다.

### WSL (Ubuntu)

```bash
cd /mnt/c/Users/gkwns/Desktop/WheelLeg/python
./.venv-wsl/bin/python main.py
```

WSLg 가 있으면(`echo $DISPLAY` 가 `:0`) matplotlib 창이 그대로 뜹니다.
SSH 접속처럼 `DISPLAY` 가 없는 환경이면 자동으로 `compare.png`, `error.png` 로 저장합니다.

### venv 재생성

```powershell
python -m venv .venv                                          # Windows
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

```bash
python3 -m venv .venv-wsl                                     # WSL / Linux
./.venv-wsl/bin/python -m pip install -r requirements.txt
```

`.venv` 폴더는 절대경로가 박혀 있어 복사/이동으로는 옮겨지지 않습니다.

## 파일 대응

| MATLAB | Python | 내용 |
|---|---|---|
| `tilde.m` | `tilde.py` | 스큐대칭 행렬 |
| `ang2mat.m` | `ang2mat.py` | RecurDyn REULER (Body 313) → 회전행렬 |
| `mat2ep.m` | `mat2ep.py` | 회전행렬 → 오일러 파라미터 (Shepperd) |
| `ep2mat.m` | `ep2mat.py` | 그 역변환 |
| `step5.m` | `step5.py` | ADAMS/RecurDyn 5차 다항 step |
| `dYdt.m` | `dYdt.py` | 상태 Y → 미분 Yp (재귀 정식화) |
| `main.m` | `main.py` | 파라미터 + RK4 + RecurDyn 비교/플롯 |
| (MATLAB 내장) | `util.py` | `col()` = `[a;b;c]`, `rms()` |

## MATLAB 과 달라지는 부분

문법상 불가피한 것만 정리했습니다.

- **인덱스**: MATLAB 1-base → Python 0-base. `Y(1:3)` → `Y[0:3]` 처럼 한 칸씩 당겨지며,
  해당 줄에 원래 MATLAB 인덱스를 주석으로 남겨두었습니다.
- **`prm` 구조체**: `types.SimpleNamespace` 로 대체. `prm.C00 = ...` 접근 방식은 동일합니다.
- **`prm.free`**: 물리적 자유도 번호라 MATLAB 과 같은 **1-base** 로 두고,
  `dYdt.py` 안에서만 0-base 로 변환합니다.
- **행렬 연산**: `*` → `@`, `'` → `.T`, `M\Q` → `np.linalg.solve(M, Q)`
- **행렬 조립**: `[A, B; C, D]` → `np.block([[A, B], [C, D]])`
- **두번째 출력**: MATLAB `[Yp, out] = dYdt(...)` → `dYdt(..., full=True)` 가 `(Yp, out)` 반환
- **플롯**: `DISPLAY` 가 없으면 PNG 저장으로 자동 전환 (MATLAB 에는 없는 분기)
- 열벡터는 shape `(n, 1)` 유지 — `(n,)` 로 두면 브로드캐스팅이 MATLAB 과 달라집니다.

## 검증

동일 초기조건에서 MATLAB 과 Python 의 전체 궤적(상태 15 + 미분 15, 2001 스텝)을 비교한 결과
**최대 상대오차 2.05e-14** (dq1(0)=0), **1.54e-14** (dq1(0)=3) 로 기계 정밀도 수준입니다.

Windows(Python 3.14.5) 와 WSL(Python 3.12.3) 실행 결과는 유효숫자 17 자리까지 완전히 동일합니다.
