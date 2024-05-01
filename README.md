
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

이 프로젝트에 기여하고 싶다면, Pull Request나 Issue를 통해 커뮤니티에 참여할 수 있습니다. 모든 기여자는 코드를 보다 나은 방향으로 발전시키는 데 도움을 줄 수 있습니다.
