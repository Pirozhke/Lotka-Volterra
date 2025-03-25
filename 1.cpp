#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>

// Структуры для параметров и состояния
struct Params { double alpha, beta, gamma, delta; };
struct State { double x, y; };

// Точка данных: время (t) и координаты (x, y)
struct DataPoint {
    double t, x, y;
    DataPoint(double _t, const State& s) : t(_t), x(s.x), y(s.y) {}
};

// Функция вычисления производных
State lotka_volterra_derivatives(const State& s, const Params& p) {
    return {
        p.alpha * s.x - p.beta * s.x * s.y,
        p.delta * s.x * s.y - p.gamma * s.y
    };
}

// Метод Рунге-Кутты 4-го порядка
State runge_kutta4(const State& s, const Params& p, double dt) {
    State k1 = lotka_volterra_derivatives(s, p);
    State k2 = lotka_volterra_derivatives({s.x + 0.5*dt*k1.x, s.y + 0.5*dt*k1.y}, p);
    State k3 = lotka_volterra_derivatives({s.x + 0.5*dt*k2.x, s.y + 0.5*dt*k2.y}, p);
    State k4 = lotka_volterra_derivatives({s.x + dt*k3.x, s.y + dt*k3.y}, p);

    return {
        s.x + (dt/6.0) * (k1.x + 2*k2.x + 2*k3.x + k4.x),
        s.y + (dt/6.0) * (k1.y + 2*k2.y + 2*k3.y + k4.y)
    };
}

// Метод Дорманда-Принса (эталон)
State dopri8(const State& s, const Params& p, double dt) {
    // Коэффициенты таблицы Бутчера (DOPRI8)
    constexpr double a[13][12] = {
        {}, // k1
        {1.0/18.0},
        {1.0/48.0, 1.0/16.0},
        {1.0/32.0, 0.0, 3.0/32.0},
        {5.0/16.0, 0.0, -75.0/64.0, 75.0/64.0},
        {3.0/80.0, 0.0, 0.0, 3.0/16.0, 3.0/20.0},
        {29443841.0/614563906.0, 0.0, 0.0, 77736538.0/692538347.0, -28693883.0/1125000000.0, 23124283.0/1800000000.0},
        {16016141.0/946692911.0, 0.0, 0.0, 61564180.0/158732637.0, 22789713.0/633445777.0, 545815736.0/2771057229.0, -180193667.0/1043307555.0},
        {39632708.0/573591083.0, 0.0, 0.0, -433636366.0/683701615.0, -421739975.0/2616292301.0, 100302831.0/723423059.0, 790204164.0/839813087.0, 800635310.0/3783071287.0},
        {246121993.0/1340847787.0, 0.0, 0.0, -37695042795.0/15268766246.0, -309121744.0/1061227803.0, -12992083.0/490766935.0, 6005943493.0/2108947869.0, 393006217.0/1396673457.0, 123872331.0/1001029789.0},
        {-1028468189.0/846180014.0, 0.0, 0.0, 8478235783.0/508512852.0, 1311729495.0/1432422823.0, -10304129995.0/1701304382.0, -48777925059.0/3047939560.0, 15336726248.0/1032824649.0, -45442868181.0/3398467696.0, 3065993473.0/597172653.0},
        {185892177.0/718116043.0, 0.0, 0.0, -3185094517.0/667107341.0, -477755414.0/1098053517.0, -703635378.0/230739211.0, 5731566787.0/1027545527.0, 5232866602.0/850066563.0, -4093664535.0/808688257.0, 3962137247.0/1805957418.0, 65686358.0/487910083.0},
        {403863854.0/491063109.0, 0.0, 0.0, -5068492393.0/434740067.0, -411421997.0/543043805.0, 652783627.0/914296604.0, 11173962825.0/925320556.0, -13158990841.0/6184727034.0, 3936647629.0/1978049680.0, -160528059.0/685178525.0, 248638103.0/1413531060.0, 0.0}
    };

    // Весовые коэффициенты 8-го порядка
    constexpr double b[] = {
        35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, 
        -2187.0/6784.0, 11.0/84.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };

    std::vector<State> k(13);
    k[0] = lotka_volterra_derivatives(s, p);

    // Вычисление стадий k2-k13
    for (int i = 1; i < 13; ++i) {
        State sum = {0.0, 0.0};
        for (int j = 0; j < i; ++j) {
            sum.x += a[i][j] * k[j].x;
            sum.y += a[i][j] * k[j].y;
        }
        State s_new = {
            s.x + dt * sum.x,
            s.y + dt * sum.y
        };
        k[i] = lotka_volterra_derivatives(s_new, p);
    }

    // Вычисление нового состояния
    State next = s;
    for (int i = 0; i < 13; ++i) {
        next.x += dt * b[i] * k[i].x;
        next.y += dt * b[i] * k[i].y;
    }

    return next;
}

// Метод Рунге-Кутты 2-го порядка (модифицированный Эйлер)
State rk2_modified_euler(const State& s, const Params& p, double dt) {
    State k1 = lotka_volterra_derivatives(s, p);
    State k2 = lotka_volterra_derivatives({
        s.x + 0.5 * dt * k1.x, 
        s.y + 0.5 * dt * k1.y
    }, p);

    return {
        s.x + dt * k2.x,
        s.y + dt * k2.y
    };
}

// Метод Рунге-Кутты 3-го порядка (Хьюна)
State rk3_heun(const State& s, const Params& p, double dt) {
    State k1 = lotka_volterra_derivatives(s, p);
    State k2 = lotka_volterra_derivatives({
        s.x + dt/3 * k1.x, 
        s.y + dt/3 * k1.y
    }, p);
    State k3 = lotka_volterra_derivatives({
        s.x + 2*dt/3 * k2.x, 
        s.y + 2*dt/3 * k2.y
    }, p);

    return {
        s.x + dt/4 * (k1.x + 3*k3.x),
        s.y + dt/4 * (k1.y + 3*k3.y)
    };
}

// Метод Рунге-Кутты 5-го порядка (Дорманда-Принса, DOPRI5)
State rk5_dormand_prince(const State& s, const Params& p, double dt) {
    // Стадия 1
    State k1 = lotka_volterra_derivatives(s, p);
    
    // Стадия 2
    State k2 = lotka_volterra_derivatives({
        s.x + dt * (1.0/5.0) * k1.x,
        s.y + dt * (1.0/5.0) * k1.y
    }, p);

    // Стадия 3
    State k3 = lotka_volterra_derivatives({
        s.x + dt * (3.0/40.0) * k1.x + dt * (9.0/40.0) * k2.x,
        s.y + dt * (3.0/40.0) * k1.y + dt * (9.0/40.0) * k2.y
    }, p);

    // Стадия 4
    State k4 = lotka_volterra_derivatives({
        s.x + dt * (44.0/45.0) * k1.x - dt * (56.0/15.0) * k2.x + dt * (32.0/9.0) * k3.x,
        s.y + dt * (44.0/45.0) * k1.y - dt * (56.0/15.0) * k2.y + dt * (32.0/9.0) * k3.y
    }, p);

    // Стадия 5
    State k5 = lotka_volterra_derivatives({
        s.x + dt * (19372.0/6561.0) * k1.x - dt * (25360.0/2187.0) * k2.x + dt * (64448.0/6561.0) * k3.x - dt * (212.0/729.0) * k4.x,
        s.y + dt * (19372.0/6561.0) * k1.y - dt * (25360.0/2187.0) * k2.y + dt * (64448.0/6561.0) * k3.y - dt * (212.0/729.0) * k4.y
    }, p);

    // Стадия 6
    State k6 = lotka_volterra_derivatives({
        s.x + dt * (9017.0/3168.0) * k1.x - dt * (355.0/33.0) * k2.x + dt * (46732.0/5247.0) * k3.x + dt * (49.0/176.0) * k4.x - dt * (5103.0/18656.0) * k5.x,
        s.y + dt * (9017.0/3168.0) * k1.y - dt * (355.0/33.0) * k2.y + dt * (46732.0/5247.0) * k3.y + dt * (49.0/176.0) * k4.y - dt * (5103.0/18656.0) * k5.y
    }, p);

    // Весовые коэффициенты для решения 5-го порядка
    return {
        s.x + dt * (35.0/384.0 * k1.x + 0.0 * k2.x + 500.0/1113.0 * k3.x + 125.0/192.0 * k4.x - 2187.0/6784.0 * k5.x + 11.0/84.0 * k6.x),
        s.y + dt * (35.0/384.0 * k1.y + 0.0 * k2.y + 500.0/1113.0 * k3.y + 125.0/192.0 * k4.y - 2187.0/6784.0 * k5.y + 11.0/84.0 * k6.y)
    };
}

// Функция для получения данных в заданные моменты времени
std::vector<DataPoint> simulate_lotka_volterra(
    const std::vector<double>& target_times,
    const State& s0,
    const Params& p,
    double dt
) {
    std::vector<DataPoint> results;
    if (target_times.empty()) return results;

    // Сортируем временные моменты по возрастанию
    std::vector<double> sorted_times = target_times;
    std::sort(sorted_times.begin(), sorted_times.end());

    State s = s0;
    double t = 0.0;
    size_t current_target = 0;

    // Пока не обработаны все целевые моменты
    while (current_target < sorted_times.size()) {
        double t_target = sorted_times[current_target];

        // Если текущее время превысило целевое, пропускаем
        if (t > t_target) {
            current_target++;
            continue;
        }

        // Вычисляем время до следующего целевого момента
        double remaining_time = t_target - t;

        // Если оставшееся время меньше шага, делаем точный шаг до цели

        s = dopri8(s, p, remaining_time);
        // s = runge_kutta4(s, p, remaining_time);
        // s = rk2_modified_euler(s, p, remaining_time);
        // s = rk3_heun(s, p, remaining_time);
        // s = rk5_dormand_prince(s, p, remaining_time);
        if (remaining_time <= dt) {
            t = t_target;
            results.emplace_back(t, s);
            current_target++;
        } 
        // Иначе делаем полный шаг dt
        else {
            t += dt;
        }
    }

    return results;
}

// Функция для записи данных графика во временный файл (остается без изменений)
void writeDataToFile(const std::vector<DataPoint>& data, const std::string& filename) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        throw std::runtime_error("Не удалось открыть файл для записи данных: " + filename);
    }

    for (const auto& point : data) {
        outputFile << point.t << " " << point.x << " " << point.y << std::endl;
    }

    outputFile.close();
}

// Функция для построения графиков
void plotGraphs(const std::vector<DataPoint>& data,
                const std::string& title, 
                const std::string& xlabel, 
                const std::string& ylabelX, 
                const std::string& ylabelY) {
    const std::string dataFile = "temp_data.dat";
    
    try {
        writeDataToFile(data, dataFile);
    } catch (...) {
        // Обработка ошибок
        return;
    }

    // Модифицированный gnuplot-скрипт
    std::string gnuplotScript = R"(
        # Первое окно: x(t) и y(t)
        set terminal wxt 0 title "Временные графики"
        set title ")" + title + R"("
        set xlabel ")" + xlabel + R"("
        set ylabel "Популяция"
        plot ")" + dataFile + R"(" using 1:2 with lines title ")" + ylabelX + R"(" , \
             ")" + dataFile + R"(" using 1:3 with lines title ")" + ylabelY + R"("

        # Второе окно: фазовый портрет y(x)
        set terminal wxt 1 title "Фазовый портрет"
        set title "Фазовый портрет"
        set xlabel ")" + ylabelX + R"("
        set ylabel ")" + ylabelY + R"("
        plot ")" + dataFile + R"(" using 2:3 with lines title "Траектория"

        # Ожидание закрытия окон
        pause mouse close
    )";

    // Запуск gnuplot
    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "%s\n", gnuplotScript.c_str());
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    }
    
    remove(dataFile.c_str());
}

int main() {
    State s0 = {15.0, 8.0};
    Params p = {0.8, 0.035, 0.3, 0.01};
    double dt = 0.0001;

    // Задаем моменты времени, для которых нужно получить данные
    std::vector<double> target_times;

    for(double i = 0; i <= 40; i += 0.0001){
        target_times.push_back(i);
    }

    // Получаем данные
    std::vector<DataPoint> trajectory = simulate_lotka_volterra(target_times, s0, p, dt);

    // Выводим результаты
    /*for (const auto& point : trajectory) {
        std::cout << "t = " << point.t 
                  << ":\tx = " << point.x 
                  << ",\ty = " << point.y << std::endl;
    }*/

    // Визуализируем графики
    plotGraphs(trajectory, "Lotka-Volterra Model", "Time (t)", "x(t)", "y(t)");

    return 0;
}