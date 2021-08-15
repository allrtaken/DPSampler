#include "RandomBits.hpp"
#include <array>
#include <algorithm>
#include <utility>
#include <functional>
#include <iostream>
#include <cassert>

RandomBits::RandomBits(): g(gmp_randinit_default) {}

string RandomBits::binary(unsigned x, uint32_t length)
{
    uint32_t logSize = (x == 0 ? 1 : log2(x) + 1);
    string s;
    do {
        s.push_back('0' + (x & 1));
    } while (x >>= 1);
    for (uint32_t i = logSize; i < (uint32_t) length; i++) {
        s.push_back('0');
    }
    std::reverse(s.begin(), s.end());

    return s;

}

string RandomBits::GenerateRandomBits(uint32_t size)
{
    string randomBits;
    std::uniform_int_distribution<unsigned> uid {0, 2147483647U};
    uint32_t i = 0;
    while (i < size) {
        i += 31;
        randomBits += binary(uid(randomEngine), 31);
    }
    randomBits.erase(size);
    return randomBits;
}

bool RandomBits::generateWeightedRandomBit(long double posWt, long double totWt){
    assert(totWt!=0);
    if (posWt == 0) return false;
    if (posWt == totWt) return true;
    long double ub = std::nextafter(totWt,std::numeric_limits<long double>::max());;
    long double lb = std::nextafter(0.0l, 1.0l);
    std::uniform_real_distribution<long double> urd(lb,ub);
    long double num = urd(randomEngine2);
    if (num>posWt){
        return false;
    } else{
        return true;
    }
}

bool RandomBits::generateWeightedRandomBit(mpf_class& posWt, mpf_class& totWt){
    assert(totWt!=0);
    if (posWt == 0) return false;
    if (posWt == totWt) return true;
    mpf_class num = g.get_f() * totWt;
    if (num>posWt){
        return false;
    } else{
        return true;
    }
}

void RandomBits::seedEngine(std::random_device& rd) {
	
    /* Initialize PRNG with seed from random_device */
    //std::random_device rd{};
    std::array<int, 10> seedArray;
    std::generate_n(seedArray.data(), seedArray.size(), std::ref(rd));
    std::seed_seq seed(std::begin(seedArray), std::end(seedArray));
    randomEngine.seed(seed);
}

void RandomBits::seedEngine2(std::random_device& rd) {
	
    /* Initialize PRNG with seed from random_device */
    //std::random_device rd{};
    std::array<int, 10> seedArray;
    std::generate_n(seedArray.data(), seedArray.size(), std::ref(rd));
    std::seed_seq seed(std::begin(seedArray), std::end(seedArray));
    randomEngine2.seed(seed);
}

void RandomBits::seedEngines() {
    std::random_device rd{};    
    seedEngine(rd);
    seedEngine2(rd);
    g.seed(rd());
}
uint64_t RandomBits::getRandInt(std::uniform_int_distribution<uint64_t>& uid){
	return uid(randomEngine);
}

double_t RandomBits::getRandReal(std::uniform_real_distribution<double_t>& urd){
	return urd(randomEngine2);
}

long double RandomBits::getRandReal(std::uniform_real_distribution<long double>& urd){
	return urd(randomEngine2);
}
