#include "excutiontrace.hpp"
#include <mutex>

using namespace std;
using namespace mcl::bn;

uint32_t calculateBlocksNeeded(uint32_t messageLengthBits) {

    const uint32_t blockBits = 512;
    messageLengthBits += 1; 
    uint32_t totalLength = messageLengthBits + 64; 


    if (totalLength % blockBits > 0) {
        totalLength += blockBits - (totalLength % blockBits);
    }


    return totalLength / blockBits;
}

size_t countOnes(const vector<Fr>& vec) {
    size_t count = 0;
    Fr one;
    one.setStr("1"); 

    for (const auto& item : vec) {
        if (item == one) { 
            ++count;
        }
    }

    return count;
}

void convertToFr(const vector<vector<vector<seed_seq::result_type>>>& input,
                 vector<vector<vector<Fr>>>& output) {
    output.resize(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        output[i].resize(input[i].size());
        for (size_t j = 0; j < input[i].size(); ++j) {
            output[i][j].resize(input[i][j].size());
            for (size_t k = 0; k < input[i][j].size(); ++k) {
                Fr element;
                input[i][j][k]; 
                output[i][j][k] = element;
            }
        }
    }
}

void serialize(const vector<vector<vector<Fr>>>& data, const string& filename) {
    ofstream file(filename, ios::binary);
    if (!file.is_open()) {
        throw runtime_error("Failed to open file for serialization");
    }

    for (const auto& matrix : data) {
        for (const auto& row : matrix) {
            for (const auto& element : row) {
                string elementStr = element.getStr();
                size_t length = elementStr.size();
                file.write(reinterpret_cast<const char*>(&length), sizeof(length));
                file.write(elementStr.data(), length);
            }
        }
    }

    file.close();
}

void deserialize(vector<vector<vector<Fr>>>& data, const string& filename) {
    ifstream file(filename, ios::binary);
    if (!file.is_open()) {
        throw runtime_error("Failed to open file for deserialization");
    }

    data.clear();
    vector<vector<Fr>> matrix;
    vector<Fr> row;
    Fr element;

    while (!file.eof()) {
        size_t length;
        file.read(reinterpret_cast<char*>(&length), sizeof(length));
        if (file.eof()) break; // 파일의 끝에 도달했는지 확인

        string elementStr(length, '\0');
        file.read(&elementStr[0], length);

        element.setStr(elementStr);
        row.push_back(element);
    }

    file.close();
}

mutex transcriptMutex; 


void commitAndAppend(KZG& PCS, polyff elem, vector<G1>& out, string& transcript, int index) {
    G1 commit_temp = PCS.commit(elem);
    out[index] = commit_temp;
    lock_guard<mutex> lock(transcriptMutex); 
    transcript += commit_temp.getStr();
}

using namespace std;
int main(int argc, char* argv[]) {

    MCL_init();

    double time=0;

    uint32_t bitSize;
    cout << "비트 사이즈를 입력하세요: ";
    cin >> bitSize;
    double b;
 
    auto start = chrono::high_resolution_clock::now();
    auto checkpoint = start;
    string transcript ="";
    uint32_t numberofblock = calculateBlocksNeeded(bitSize);
    uint32_t domain_size = nextPowOf2_32(6856*numberofblock + 3*(1<<12));
    KZG PCS;
    PCS.init(domain_size+5);
    Ntt<Fr> ntt;
    ntt.init(domain_size);
    const Fr* ptrW = ntt.getWs();
    vector<Fr> Ws(ptrW, ptrW + domain_size);
    vector<Fr> tweedle(domain_size, Fr(0));
    tweedle.push_back(Fr(1)); 
    tweedle[0] = -Fr(1);
    polyff cyclotomic;
    cyclotomic.data = tweedle;
    polyff cyclotomicwx = cyclotomic.composeWithMonomial(Ws[1]);
    polyff numerator  = cyclotomic.syntheticDivide(Fr(-1));
    const polyff L_0 = numerator *(1/numerator.evaluate(Ws[0]));
    circuit A;
    auto end_setup = chrono::high_resolution_clock::now();
    cout << "설정 소요시간 : " << duration_cast<milliseconds>(end_setup - checkpoint).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();
    

    tie(A, b) = generateExcutiontrace(bitSize);
    time += b;
    A.excutiontrace.data.resize(domain_size,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});
    A.init_table("tablecsv_backup.csv");

    auto end_trace = chrono::high_resolution_clock::now();
    cout << "실행추적 생성 소요시간 : " << duration_cast<milliseconds>(end_trace - checkpoint).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();

    vector<Fr> selector_tmp(domain_size,0);
    polyff selector_poly_tmp;
    for(int i = 0; i < 3; i++){

        selector_tmp = A.excutiontrace.getColumn(i);
        ntt.intt(&selector_tmp[0]);
        selector_poly_tmp.data = selector_tmp;
        A.poly_selector.push_back(selector_poly_tmp);

    }

    auto position = A.permutation_pos_interior;
    for(int i = 0; i < position.size(); i++){

        position[i][0][1] -= 3;
        position[i][1][1] -= 3;

    }
    auto[S, S_ij, index_origin] = permutation(Ws, domain_size, 7, position);
    A.poly_permutation = S;

    auto end_pre = chrono::high_resolution_clock::now();
    cout << "테이블, 선택자, 순열 다항식 생성 소요시간 : " << duration_cast<milliseconds>(end_pre - checkpoint).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();

    Fr b1, b2;
    matrixFr AA(domain_size, 7);
    for(int i = 0; i < 7; i++){

        AA.replaceColumn(i, A.excutiontrace.getColumn(i+3));

    }

    AA.duplicateMultipleElements(position);


    vector<Fr> wire_tmp(domain_size,0);
    polyff wire_poly_tmp;
    for(int i = 0; i < 7; i++){

        b1.setByCSPRNG();
        b2.setByCSPRNG();

        wire_tmp = AA.getColumn(i);
        ntt.intt(&wire_tmp[0]);
        wire_poly_tmp.data = wire_tmp;
        wire_poly_tmp = polyff({b2, b1})*cyclotomic + wire_poly_tmp;
        A.poly_wire.push_back(wire_poly_tmp);

    }

    for(int i = 10; i < 16; i++){

        b1.setByCSPRNG();
        b2.setByCSPRNG();

        wire_tmp = A.excutiontrace.getColumn(i);
        ntt.intt(&wire_tmp[0]);
        wire_poly_tmp.data = wire_tmp;
        wire_poly_tmp = polyff({b2, b1})*cyclotomic + wire_poly_tmp;
        A.poly_wire.push_back(wire_poly_tmp);

    }

    AA.duplicateMultipleElements(position);
    auto end_wire = chrono::high_resolution_clock::now();
    cout << "와이어 다항식 생성 소요시간 : " << duration_cast<milliseconds>(end_wire - checkpoint).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();
    

    vector<thread> threads;
    vector<G1> out1st(A.poly_wire.size());
    G1 commit_temp;

    for (int i = 0; i < A.poly_wire.size(); ++i) {
        threads.emplace_back(commitAndAppend, ref(PCS), A.poly_wire[i], ref(out1st), ref(transcript), i);
    }

    for (auto& t : threads) {
        t.join();
    }

    Fr zeta;
    zeta.setHashOf(transcript);

    auto end_tscr = chrono::high_resolution_clock::now();
    cout << "랜덤오라클 쿼리1 소요시간 : " << duration_cast<milliseconds>(end_tscr - checkpoint).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();

    auto[f, t, tid] = lookup(A.excutiontrace, A.table, 2, 3, 5, zeta);
    auto tid_copy = tid;
    ntt.intt(&tid_copy[0]);
    polyff f_tid;
    f_tid.data = tid_copy;
    auto[fx, tx, twx, h1x, h1wx, h2x, h2wx] = plookup_copyconstraints_1(f, t);

    Fr b3, b4, b5, b6, b7;
    b3.setByCSPRNG();
    b4.setByCSPRNG();
    b5.setByCSPRNG();
    b6.setByCSPRNG();
    b7.setByCSPRNG();

    A.fx = polyff({b2, b1})*cyclotomic + fx;
    A.tx = tx;
    A.twx = twx;
    A.h1x = polyff({b5, b4, b3})*cyclotomic + h1x;
    A.h1wx = polyff({b5, b4, b3}).composeWithMonomial(Ws[1])*cyclotomicwx + h1wx;
    A.h2x = polyff({b7, b6})*cyclotomic + h2x;
    A.h2wx = polyff({b7, b6}).composeWithMonomial(Ws[1])*cyclotomicwx + h2wx;

    auto end_plookup = chrono::high_resolution_clock::now();
    cout << "plookup 소요시간 : " << duration_cast<milliseconds>(end_plookup - checkpoint).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();

    vector<G1> out2nd(3); 
    vector<polyff> transcript2_poly = {fx, h1x, h2x}; 


    vector<thread> commitThreads;
    for (int i = 0; i < transcript2_poly.size(); ++i) {
        commitThreads.emplace_back(commitAndAppend, ref(PCS), transcript2_poly[i], ref(out2nd), ref(transcript), i);
    }


    for (auto& t : commitThreads) {
        t.join();
    }


    for (const auto& commit_result : out2nd) {
        transcript += commit_result.getStr();
    }

    Fr beta1, gamma1, delta1, epsillon1;
    beta1.setHashOf(transcript+"1");
    gamma1.setHashOf(transcript+"2");
    delta1.setHashOf(transcript+"3");
    epsillon1.setHashOf(transcript+"4");

    auto end_tscr1 = chrono::high_resolution_clock::now();
    cout << "랜덤오라클 쿼리2 소요시간 : " << duration_cast<milliseconds>(end_tscr1 - checkpoint).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();

    auto result = Plookup_poly_z(f, t, beta1, gamma1);
    polyff zx = get<0>(result); 
    polyff zwx = get<1>(result);


    polyff accx = copyconstraints(AA, S, S_ij, index_origin, delta1, epsillon1, position);
    polyff accwx = accx.composeWithMonomial(Ws[1]);

    Fr b8, b9, b10;
    b8.setByCSPRNG();
    b9.setByCSPRNG();
    b10.setByCSPRNG();

    A.zx = polyff({b10, b9, b8})*cyclotomic + zx;
    A.zwx = polyff({b10, b9, b8}).composeWithMonomial(Ws[1])*cyclotomicwx + zwx;

    
    polyff one = {1};
    polyff beta = {beta1};
    polyff gamma = {gamma1};

    polyff tz_z2 = A.zwx*((gamma*(one+beta)+A.h1x+beta*A.h1wx)*(gamma*(one+beta)+A.h2x+beta*A.h2wx))
    -A.zx*((one+beta)*(gamma+A.fx)*(gamma*(one+beta)+A.tx+beta*A.twx));
    polyff tz_L0_z2 = (zx-1)*L_0;

    auto end_z2 = chrono::high_resolution_clock::now();
    cout << "z_2(x)생성 소요시간 : " << duration_cast<milliseconds>(end_z2 - checkpoint).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();    

    Fr b11, b12, b13;
    b11.setByCSPRNG();
    b12.setByCSPRNG();
    b13.setByCSPRNG();

    A.accx = polyff({b13, b12, b11})*cyclotomic + accx;
    A.accwx = polyff({b13, b12, b11}).composeWithMonomial(Ws[1])*cyclotomicwx + accwx;


    polyff delta = {delta1};
    polyff epsillon = {epsillon1};
    polyff zlhs;
    zlhs = accwx
    *((A.poly_wire[0]+delta*S[0]+epsillon)
    *(A.poly_wire[1]+delta*S[1]+epsillon)
    *(A.poly_wire[2]+delta*S[2]+epsillon)
    *(A.poly_wire[3]+delta*S[3]+epsillon)
    *(A.poly_wire[4]+delta*S[4]+epsillon)
    *(A.poly_wire[5]+delta*S[5]+epsillon)
    *(A.poly_wire[6]+delta*S[6]+epsillon));



    polyff zrhs;
    zrhs = accx
    *((A.poly_wire[0]+delta*polyff({0,1})+epsillon)
    *(A.poly_wire[1]+delta*polyff({0,2})+epsillon)
    *(A.poly_wire[2]+delta*polyff({0,3})+epsillon)
    *(A.poly_wire[3]+delta*polyff({0,4})+epsillon)
    *(A.poly_wire[4]+delta*polyff({0,5})+epsillon)
    *(A.poly_wire[5]+delta*polyff({0,6})+epsillon)
    *(A.poly_wire[6]+delta*polyff({0,7})+epsillon));

    polyff tz_z1 = zlhs.PolyCondense() - zrhs.PolyCondense();
    polyff tz_L0_z1 = (accx-1)*L_0;
    tz_z1 = tz_z1.PolyCondense();

    auto end_z1 = chrono::high_resolution_clock::now();
    cout << "z_1(x)생성 소요시간 : " << duration_cast<milliseconds>(end_z1 - checkpoint).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();

    vector<G1> out3rd(2); 
    vector<polyff> commitPolys = {zx, accx};


    vector<thread> commitThreads2;


    for (int i = 0; i < commitPolys.size(); ++i) {
        commitThreads2.emplace_back(commitAndAppend, ref(PCS), commitPolys[i], ref(out3rd), ref(transcript), i);
    }


    for (auto& t : commitThreads2) {
        t.join();
    }


    for (const auto& commit_result : out3rd) {
        transcript += commit_result.getStr();
    }

    Fr alpha1;
    alpha1.setHashOf(transcript);
    
    polyff alpha = {alpha1};
    polyff alpha_2 = {alpha1*alpha1};
    polyff alpha_3 = {alpha1*alpha1*alpha1};
    polyff alpha_4 = {alpha1*alpha1*alpha1*alpha1};
    polyff alpha_5 = {alpha1*alpha1*alpha1*alpha1*alpha1};
        

    auto end_tscr2 = chrono::high_resolution_clock::now();
    cout << "랜덤오라클 쿼리3 소요시간 : " << duration_cast<milliseconds>(end_tscr2 - checkpoint).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();

    polyff tz1;

    tz1 = A.poly_selector[0]
    *((A.poly_wire[0]*A.poly_wire[7])
    +(A.poly_wire[1]*A.poly_wire[8])
    +(A.poly_wire[2]*A.poly_wire[9])
    +(A.poly_wire[3]*A.poly_wire[10])
    +(A.poly_wire[4]*A.poly_wire[11])
    +(A.poly_wire[5]*A.poly_wire[12])-A.poly_wire[6]);

    auto end_tz1 = chrono::high_resolution_clock::now();
    cout << "t(x)z_H(X) 제약조건1 생성 소요시간 : " << duration_cast<milliseconds>(end_tz1 - checkpoint).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();    

    polyff tz2;
    tz2 = A.poly_selector[1]
    *(A.poly_wire[0]+A.poly_wire[1]+A.poly_wire[2]+A.poly_wire[3]+A.poly_wire[4]+A.poly_wire[5]-A.poly_wire[6]);

    auto end_tz2 = chrono::high_resolution_clock::now();
    cout << "t(x)z_H(X) 제약조건2 생성 소요시간 : " << duration_cast<milliseconds>(end_tz2 - checkpoint).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();

    polyff zetap = {zeta};
    polyff zetap2 = {zeta*zeta};
    polyff zetap3 = {zeta*zeta*zeta};
    polyff zetap4 = {zeta*zeta*zeta*zeta};
    A.fx_t  =   polyff({b2, b1})*cyclotomic + A.poly_wire[0]+zetap*A.poly_wire[1]+zetap2*A.poly_wire[2]+zetap3*A.poly_wire[3]
    +zetap4*f_tid;
    polyff tz3 = A.poly_selector[2]
    *(A.poly_wire[0]+zetap*A.poly_wire[1]+zetap2*A.poly_wire[2]+zetap3*A.poly_wire[3]
    +zetap4*f_tid - A.fx_t);

    auto end_tz3 = chrono::high_resolution_clock::now();
    cout << "t(x)z_H(X) lookup 제약조건 생성 소요시간 : " << duration_cast<milliseconds>(end_tz3 - checkpoint).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();

    A.tz = tz1 + tz2 + alpha*tz_z1 + alpha_2*tz_L0_z1 + alpha_3*tz3 + alpha_4*tz_z2 + alpha_5*tz_L0_z2;

    auto end = chrono::high_resolution_clock::now();
    cout << "총 소요시간(설정 제외) : " << duration_cast<milliseconds>(end - end_setup).count() << " 밀리초\n";
    checkpoint = chrono::high_resolution_clock::now();

    return 0;
}
