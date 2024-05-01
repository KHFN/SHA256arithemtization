
### 필수 조건

프로젝트를 빌드하고 실행하기 전에 다음 소프트웨어와 라이브러리가 시스템에 설치되어 있어야 합니다:

- GCC 컴파일러 (C++11 이상 지원)
- CMake (선택적, 몇몇 빌드 스크립트에 필요)
- Git
- Python 3 (빌드 스크립트 실행에 사용)
- [herumi/mcl](https://github.com/herumi/mcl) 라이브러리: 고성능 암호 연산 라이브러리

### 설치

1. **리포지토리 복제**

   GitHub에서 프로젝트를 복제합니다:

   ```bash
   git clone https://github.com/KHFN/SHA256arithemtization.git
   cd SHA256arithemtization
   ```

2. **서브모듈 초기화 및 업데이트**

   프로젝트에는 `herumi/mcl`과 같은 필수 서브모듈이 포함되어 있습니다. 다음 명령으로 초기화 및 업데이트를 진행합니다:

   ```bash
   git submodule update --init --recursive
   ```

3. **빌드**

   빌드를 위해 기본 디렉토리에서 다음 명령을 실행합니다:

   ```bash
   cd mcl
   make -j4
   ```

   이 명령은 멀티 코어 시스템에서 병렬로 컴파일을 수행하여 시간을 절약할 수 있습니다.

### 실행

빌드가 완료된 후, `arithmetization` 폴더 내의 `main.cpp` 파일을 컴파일하여 실행 파일을 생성합니다:

```bash
cd ../arithmetization
g++ -o main_fast main.cpp -I../mcl/include/ -L../mcl/lib -lmcl -lpthread -Wl,-rpath,../mcl/lib -Ofast
./main_fast
```

## Android 스마트폰에서 실행해보기

Termux는 Android 기기에서 리눅스 환경을 제공하는 앱입니다. 이를 통해 GCC, Python 등 리눅스 기반 도구를 사용할 수 있습니다.

### 필수 조건

- Android 기기
- Termux 앱 (F-Droid에서 다운로드 가능)

### 설치 과정

1. **Termux 설치**

   Google Play Store 또는 F-Droid에서 Termux 앱을 설치합니다.

2. **필수 패키지 업데이트 및 설치**

   Termux를 열고 다음 명령으로 시스템을 업데이트하고 필요한 도구를 설치합니다:

   ```bash
   cd ..
   pkg update && pkg upgrade
   pkg install git python clang make cmake -y
   ```

3. **리포지토리 복제**

   다음 명령으로 프로젝트 리포지토리를 복제합니다:

   ```bash
   git clone https://github.com/KHFN/SHA256arithemtization.git
   cd SHA256arithemtization
   ```

4. **서브모듈 초기화 및 업데이트**

   프로젝트의 서브모듈을 초기화하고 업데이트합니다:

   ```bash
   git submodule update --init --recursive
   ```

5. **mcl 라이브러리 빌드**

   mcl 라이브러리를 빌드하기 전에 `CXX` 환경 변수를 `clang++`로 설정하여 `clang` 컴파일러를 사용합니다:

   ```bash
   cd mcl
   mkdir build
   cd build
   cmake .. -DCMAKE_CXX_COMPILER=clang++
   make
   cp -r lib ../
   ```

### 프로젝트 실행

빌드가 완료된 후, `arithmetization` 디렉토리로 돌아가 `main.cpp` 파일을 컴파일하고 실행 파일을 생성합니다:

```bash
cd ../../arithmetization
clang++ -o main_fast main.cpp -I../mcl/include/ -L../mcl/lib -lmcl -lpthread -Wl,-rpath,../mcl/lib -Ofast
./main_fast
```
