#ifndef FTLIB_H
#define FTLIB_H

#include <vector>
#include <complex>

typedef unsigned char uchar;

constexpr const double pi = 3.141592653589793238463;

static const uchar rtable[256] = {
	0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
	0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
	0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
	0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
	0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
	0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
	0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
	0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
	0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
	0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
	0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
	0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
	0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
	0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
	0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
	0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
};


// @TODO: delete this function
#include <iostream>
void print_bits(size_t num)
{
	int bit(sizeof(size_t) * CHAR_BIT - 1), counter(0);
	for (; bit >= 0; bit--)
		std::cout << (size_t)((num >> bit) & 1);
}

template <typename UINT>
UINT reverse(UINT unum, uchar bits = sizeof(UINT)*CHAR_BIT)  
{
	static_assert(std::is_unsigned<UINT>::value, 
		"reverse only unsigned type");

	UINT result;
	uchar *input = reinterpret_cast<uchar*>(&unum);
	uchar *output = reinterpret_cast<uchar*>(&result);
	for (size_t i = 0; i < sizeof(UINT); i++)
		output[i] = rtable[input[sizeof(UINT) - i - 1]];
	return result >> (sizeof(UINT)*CHAR_BIT - bits);
}


template < typename SignalType >
std::vector<std::complex<SignalType>> dft_formula(const std::vector<SignalType>& signal)
{
	const auto N = signal.size();

	std::vector<std::complex<SignalType>> result(N);

	const SignalType k = -2.0 * pi / N;

	for (size_t n = 0; n < N; n++) {
		for (size_t m = 0; m < N; m++) {
			std::complex<SignalType> X {
				cos(k * n * m),
				sin(k * n * m)
			};
			result[n] += signal[m] * X;
		}
	}

	return result;
}

template < typename SignalType >
std::vector<std::complex<SignalType>> idft_formula(const std::vector<std::complex<SignalType>> &spectr)
{
	const auto N = spectr.size();
	const SignalType k = 2.0 * pi / N;
	std::vector<std::complex<SignalType>> result(N);

	for (size_t n = 0; n < N; n++) {
		for (size_t m = 0; m < N; m++) {
			std::complex<SignalType> X{
				cos(k * n * m),
				sin(k * n * m)
			};
			result[n] += spectr[m] * X;
		}
		result[n] /= N;
	}


	return result;
}

std::vector<std::complex<float>> fdft_table(const std::vector<float> &signal)
{
	const size_t N = signal.size();
	const uchar bits = log2(N);

	std::vector<std::complex<float>> result(N);

	for (size_t i = 0; i < N; i++) {
		const size_t idx = reverse<size_t>(i, bits);
		result[i] = std::complex<float>(signal[idx]);
	}

	for (size_t len_sub_array = 2; len_sub_array <= N;) {
		const float k = -2.0 * pi / len_sub_array;
		for (size_t sub_array_idx = 0; sub_array_idx < N;) {
			for (size_t i = sub_array_idx; i < sub_array_idx + len_sub_array / 2; i++) {
				std::complex<float> W { cos(k * i), sin(k * i) };
				auto tmp = result[i] + W * result[i + len_sub_array / 2];
				result[i + len_sub_array / 2] = result[i] - W * result[i + len_sub_array / 2];
				result[i] = tmp;
			}
			sub_array_idx += len_sub_array;
		}
		len_sub_array <<= 1;
	}
	
	return result;
}


bool fdft_table(const float *signal, std::complex<float> *result, const size_t N)
{
	const uchar bits = log2(N);

	for (size_t i = 0; i < N; i++) {
		const size_t idx = reverse<size_t>(i, bits);
		result[i] = std::complex<float>(signal[idx]);
	}

	for (size_t len_sub_array = 2; len_sub_array <= N;) {
		const float k = -2.0 * pi / len_sub_array;
		for (size_t sub_array_idx = 0; sub_array_idx < N;) {
			for (size_t i = sub_array_idx; i < sub_array_idx + len_sub_array / 2; i++) {
				std::complex<float> W{ cos(k * i), sin(k * i) };
				auto tmp = result[i] + W * result[i + len_sub_array / 2];
				result[i + len_sub_array / 2] = result[i] - W * result[i + len_sub_array / 2];
				result[i] = tmp;
			}
			sub_array_idx += len_sub_array;
		}
		len_sub_array <<= 1;
	}

	return true;
}




template < typename SignalType >
std::vector<std::complex<SignalType>> fdft_recursive(const std::vector<SignalType> &signal)
{
	const size_t N = signal.size();

	std::vector<std::complex<SignalType>> result(N);

	if (N <= 2) {
		if (N == 2) {
			result[0] = std::complex<SignalType>(signal[0] + signal[1]);
			result[1] = std::complex<SignalType>(signal[0] - signal[1]); 
			return result;
		}
		if (N == 1) result[0] = std::complex<SignalType>(signal[0]);
		return result;
	}

	std::vector < SignalType > even(N / 2);
	std::vector< SignalType > odd(N / 2);

	for (size_t i = 0; i < N / 2; i++) {
		even[i] = signal[2*i];
		odd[i] = signal[2*i + 1];
	}

	auto evenRes = fdft(even);
	auto oddRes = fdft(odd);

	const SignalType k = -2.0 * pi / N;

	for (size_t n = 0; n < N / 2; n++) {
		std::complex<SignalType> a(cos(k * n), sin(k * n));
		result[n] = evenRes[n] + a*oddRes[n];
		result[n + N / 2] = evenRes[n] - a*oddRes[n];
	}

	return result;
}

#endif FTLIB_H 
