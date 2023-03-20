#pragma once
//고려대학교 김원이 만든 Elgamal 암호입니다.
//Elgamal의 KeyGen/Encryption/Decryption을 GMP library를 이용해 구현한 프로그램입니다.

#ifndef __ELGAMAL_HPP__
#define __ELGAMAL_HPP__

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "gmp.h"
#include "gmpxx.h"

std::vector<mpz_class> sieve_of_Eratosthenes(mpz_class n);
std::vector<mpz_class> primes_lowerthan(mpz_class n);

mpz_class Pollard_Rho(mpz_class n); //integer factorizing
mpz_class primitive_root(mpz_class p); //primitive root generator

bool subsafe_check(mpz_class n);

std::vector<mpz_class> keyGen(std::vector<mpz_class>& v);
std::vector<mpz_class> Encrypt(mpz_class m, mpz_class g, mpz_class h, mpz_class p);
mpz_class Decrypt(mpz_class c1, mpz_class c2, mpz_class d, mpz_class p);

#endif

/*
References

1. Zarni Sann, Thi Thi Soe, Khaing Myat Nwe, Comparison of Public Key Cryptography in Different Security Level, International Journal of Recent Development in Engineering and Technology Volume 8, Issue 12, December 2019
2. Behrouz Forouzan, "Cryptography and Network Security", 3rd Edition
3. Serge Vaudenay, "A Classical Introduction to Cryptography", 2006
4. Kenneth H. Rosen, "Elementary Number Theory and Its Application", 6th Edition
5. Henri Cohen, "A Course in Computational Algebraic Number Theory", 2nd Edition
6. GNU MP 6.2.1 Manual
7. Michael J. Wiener, Safe Prime Generation with a Combined Sieve, Cryptographic Clarity, 20 Hennepin St., Nepean, Ontario, Canada K2J 3Z4

Not Used
8. AkinTunde Ayodele, Here's How Quadratic Sieve Factorization Works, Nerd For Tech, Feburary 6, 2022
9. A. Kruppa, Comparison of Block-Lanczos and Wiedemann for Solving Linear Systems in Large Factorizations, Workshop on Computational Number Theory 2011
10. Sonia Belaid, Sylvain Lachartre, Solving Sparse Systems with the Block Wiedemann Algorithm, Seminaire SALSA 2011, 8 July
11. Joachim von zur Gathen, Igor E. Shparlinski, Generating Safe Primes, Journal of Mathematical Cryptology Volume 7, Issue 4, 2013
*/