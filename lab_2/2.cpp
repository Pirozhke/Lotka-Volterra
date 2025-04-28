#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Структура для хранения числа в формате double-double
typedef struct {
    double hi;  // Старшие биты числа
    double lo;  // Младшие биты числа (|lo| ≤ 0.5 * ulp(hi))
} ddouble;

//--------------------------------------------------
// Арифметические операции с двойной точностью
//--------------------------------------------------

// Алгоритм Деккера для точного сложения
ddouble dd_add(ddouble a, ddouble b) {
    // Вычисление старшей части суммы и погрешности
    double s1 = a.hi + b.hi;
    double s1_err = s1 - a.hi;
    
    // Компенсация младших разрядов
    double s2 = (b.hi - s1_err) + (a.hi - (s1 - s1_err)) + a.lo + b.lo;
    
    // Объединение результатов
    return (ddouble){s1 + s2, s2 - ((s1 + s2) - s1)};
}

// Точное умножение с использованием FMA
ddouble dd_mul(ddouble a, ddouble b) {
    // Разложение произведения на старшую и младшую части
    double p1 = a.hi * b.hi;
    double p2 = fma(a.hi, b.hi, -p1);  // Точное вычисление ошибки
    
    // Учет вклада младших частей
    p2 = fma(a.hi, b.lo, p2);
    p2 = fma(a.lo, b.hi, p2);
    p2 = fma(a.lo, b.lo, p2);
    
    // Объединение результатов
    return (ddouble){p1 + p2, p2 - ((p1 + p2) - p1)};
}

// Константа ln(2) в формате double-double
static const ddouble dd_ln2 = {
    0.693147180559945309417232121458176568,  // Старшие биты
    2.3190468138462995584177715810474248e-17 // Младшие биты
};

//--------------------------------------------------
// Вычисление натурального логарифма
//--------------------------------------------------
ddouble dd_ln(ddouble x) {
    assert(x.hi > 0 && "Ошибка: логарифм отрицательного числа");
    
    // Нормализация аргумента: x = m * 2^exp, где m ∈ [0.5, 1)
    int exp;
    double m_hi = frexp(x.hi, &exp);
    ddouble m = {m_hi, ldexp(x.lo, -exp)};
    
    // Разложение в ряд Тейлора для ln(1 + t)
    ddouble t = dd_add(m, (ddouble){-1.0, 0.0});
    ddouble sum = t;
    ddouble term = t;

    // Итеративное вычисление суммы ряда
    for (int k = 2; k <= 100; k++) {
        term = dd_mul(term, t);
        ddouble coeff = {(k % 2 == 0) ? -1.0/k : 1.0/k, 0.0};
        sum = dd_add(sum, dd_mul(term, coeff));
    }

    // Учет экспоненты и добавление ln(2^exp)
    ddouble exp_term = dd_mul((ddouble){(double)exp, 0.0}, dd_ln2);
    return dd_add(exp_term, sum);
}

//--------------------------------------------------
// Вычисление экспоненты
//--------------------------------------------------
ddouble dd_exp(ddouble x) {
    // Редукция аргумента: x = k*ln2 + r, где |r| ≤ ln2/2
    double k_d = rint((x.hi + x.lo) / dd_ln2.hi);
    ddouble k = {k_d, 0.0};
    ddouble kln2 = dd_mul(k, dd_ln2);
    ddouble r = dd_add(x, dd_mul(kln2, (ddouble){-1.0, 0.0}));

    // Разложение в ряд Тейлора для e^r
    ddouble sum = {1.0, 0.0};
    ddouble term = {1.0, 0.0};
    
    for (int n = 1; n <= 100; n++) {
        term = dd_mul(term, r);
        term = dd_mul(term, (ddouble){1.0/n, 0.0});
        sum = dd_add(sum, term);
    }

    // Восстановление результата: e^x = 2^k * e^r
    double two_k = pow(2.0, k_d);
    return dd_mul(sum, (ddouble){two_k, 0.0});
}

//--------------------------------------------------
// Возведение в степень через экспоненту и логарифм
//--------------------------------------------------
ddouble dd_pow(ddouble x, ddouble y) {
    ddouble ln_x = dd_ln(x);
    ddouble exponent = dd_mul(y, ln_x);
    return dd_exp(exponent);
}

//--------------------------------------------------
// Визуализация различий с эталонным значением
//--------------------------------------------------
void print_diff(const char* value, const char* reference) {
    for (int i = 0; i < 40 && reference[i]; i++) {
        if (value[i] == reference[i]) {
            printf("%c", value[i]);
        } else {
            printf("\033[31m%c\033[0m", value[i]); // Красный цвет для различий
        }
    }
    printf("\n");
}

// Форматированный вывод сравнения результатов
void print_comparison(const char* label, double x, double y, const char* wolfram) {
    ddouble dd = dd_pow((ddouble){x, 0.0}, (ddouble){y, 0.0});
    char d_str[128], dd_str[128];
    
    // Преобразование результатов в строки
    snprintf(d_str, sizeof(d_str), "%.40f", pow(x, y));
    snprintf(dd_str, sizeof(dd_str), "%.40f", dd.hi + dd.lo);
    
    // Вывод сравнения
    printf("\n\033[1m%s\033[0m\n", label);
    printf("Double:     "); print_diff(d_str, wolfram);
    printf("DDouble:    "); print_diff(dd_str, wolfram);
    printf("Wolfram:    \033[36m%s\033[0m\n", wolfram);
}

//--------------------------------------------------
// Тестовые примеры для проверки точности
//--------------------------------------------------
int main() {
    print_comparison("Пример 1: 2^3", 2, 3, "8.0000000000000000000000000000000000000000");
    print_comparison("Пример 2: 3^2", 3, 2, "9.0000000000000000000000000000000000000000");
    print_comparison("Пример 3: 2.5^1.5", 2.5, 1.5, "3.9528470752104740643519638617963167323809");
    print_comparison("Пример 4: π^π", M_PI, M_PI, "36.46215960720791177099082602269212366636");
    print_comparison("Пример 5: 1.7^4.3", 1.7, 4.3, "9.793329391416697851965838346547452611583");
    
    return 0;
}