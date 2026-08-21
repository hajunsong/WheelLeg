# python — matlab/ 의 1:1 이식

MATLAB 코드와 파일명·함수명·변수명을 그대로 유지해 옮긴 것입니다.
Windows / WSL / 네이티브 리눅스에서 모두 같은 결과가 나옵니다.

## 실행

플랫폼별로 venv 를 따로 둡니다 (`Scripts/*.exe` 와 `bin/` 구조가 호환되지 않습니다).

| 플랫폼 | venv | 셋업 |
|---|---|---|
| Windows | `.venv` | `python -m venv .venv` + `pip install -r requirements.txt` |
| Linux / WSL / macOS | `.venv-linux` | `./setup.sh` |

### Linux / WSL / macOS

```bash
cd python
./setup.sh                          # venv 생성 + 의존성 설치
./.venv-linux/bin/python main.py
./.venv-linux/bin/python lqr_control.py
```

`PYTHON=python3.12 ./setup.sh` 처럼 인터프리터를 지정할 수도 있습니다.

`venv` 모듈이 없다는 오류가 나면 배포판 패키지를 먼저 설치하세요.

```bash
sudo apt install python3-venv        # Debian / Ubuntu
sudo dnf install python3-virtualenv  # Fedora / RHEL
```

**시각화**: `DISPLAY` 나 `WAYLAND_DISPLAY` 가 있으면 창이 뜹니다 (WSL 은 WSLg 로 자동).
없으면(SSH, 컨테이너, CI 등) PNG로 저장합니다. `lqr_control.py`는 실행 환경과
관계없이 `lqr_animation.gif`도 생성합니다.

### Windows (PowerShell)

```powershell
cd C:\Users\gkwns\Desktop\WheelLeg\python
.\.venv\Scripts\python.exe main.py
.\.venv\Scripts\python.exe lqr_control.py
```

처음이라면 venv 부터:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

활성화해서 쓰려면 `.\.venv\Scripts\Activate.ps1` 실행 후 `python main.py`, 나올 때 `deactivate`.
VSCode 에서는 `Ctrl+Shift+P` → **Python: Select Interpreter** 로
`.\python\.venv\Scripts\python.exe` 를 지정하면 실행 버튼이 그대로 동작합니다.

> venv 폴더는 내부에 절대경로가 박혀 있어 복사·이동으로는 옮겨지지 않습니다. 항상 새로 만드세요.

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
| — | `model.py` | SI 모델 파라미터 + 축약/전체 상태 변환 |
| — | `lqr_control.py` | 수치 선형화 + CARE 기반 LQR + 비선형 폐루프 시뮬레이션 |
| (MATLAB 내장) | `util.py` | `col()` = `[a;b;c]`, `rms()` |
| — | `setup.sh` | POSIX 환경 셋업 (Python 전용) |

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
- **제어 입력**: `dYdt(..., Fy=u)`로 시간 함수 대신 외력 `u`[N]를 직접 입력
- **플롯**: 디스플레이가 없으면 PNG 저장으로 자동 전환 (MATLAB 에는 없는 분기)
- 열벡터는 shape `(n, 1)` 유지 — `(n,)` 로 두면 브로드캐스팅이 MATLAB 과 달라집니다.

## 이식성

- **경로**: 기준 데이터 `rec_data.csv` 를 `Path(__file__)` 기준으로 찾으므로
  작업 디렉터리와 무관하게 실행됩니다.
- **대소문자**: 네이티브 리눅스(ext4)는 대소문자를 구분합니다. `/mnt/c` 나 NTFS 에서는
  구분하지 않아 파일명 불일치가 드러나지 않으므로, 검증은 ext4 에 clone 해서 수행했습니다.
- **줄바꿈**: 저장소 루트 `.gitattributes` 로 `.py` 는 LF 고정입니다.

## 검증

- MATLAB ↔ Python : 동일 초기조건에서 전체 궤적(상태 15 + 미분 15, 2001 스텝) 비교 결과
  **최대 상대오차 2.05e-14**
- Windows (Python 3.14.5) ↔ Linux (Python 3.12.3) : 유효숫자 **17 자리까지 완전 일치**
- 네이티브 리눅스 검증 : GitHub 에서 ext4 로 clone → `./setup.sh` → 실행까지 확인
