#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>


// Входной файл с координатами линий
#define INPUT_NAME "data\\test_files\\3-lines.txt"
// Название выходного для тестирования (для проверки промежуточных результатов)
#define TESTING_OUTPUT_NAME "data\\test_output.txt"
// Название выходного файла (файл, куда алгоритм выдаст точки)
#define OUTPUT_NAME "data\\output.txt"
// Минимальные и максимальные координаты линий
#define XMIN -2000
#define YMIN -2000
#define XMAX 2000
#define YMAX 2000
// Левый верхний угол сетки
#define PMIN_X 990
#define PMIN_Y 990
// Кол-во боксов в сетке
#define NBOXES_X 1010
#define NBOXES_Y 1010
// Размер боксов в сетке
// Выставлять одинаковыми, ибо функция расчитывания расстояния оптимизирована
//  благодаря предположению, что ys == xs.
// Если хотите использовать разные, напишите свою функцию расчитывания расстояния.
#define XS 2
#define YS 2
// Множитель для вычисления кол-ва top боксов (TOP_1 * N)
#define TOP_1 4 
// Максимальное время должно быть больше,
// чем расстояние до границы деленное длину линии
// Условие:
//      |MAX_TIME * (x2 - x1)| > |XMAX - XMIN|  AND
//      |MAX_TIME * (y2 - y1)| > |YMAX - YMIN| 
#define MAX_TIME 10000000.0
// Минимальное расстояние(евклидово) между кластерами
#define D_EU 50
// D = (D_EU / box_size)^2, где D_EU - обычное евклидовое расстояние 
const double D = std::pow(double(D_EU) / XS, 2);   


/**
 * @brief Создает динамическую матрицу
 * 
 * @tparam T 
 * @param m 
 * @param n 
 * @param val 
 * @return T** 
 */
template<class T>
T** create_matrix(const short m, const short n, const T& val) {
    T** mat = new T*[m];
    for (short i = 0; i < m; i++) {
        mat[i] = new T[n];
        for (short j = 0; j < n; j++) {
            mat[i][j] = val;
        }
    }
    return mat;
}


/**
 * @brief Удаляет динамическую матрицу
 * 
 * @tparam T 
 * @param mat 
 * @param m 
 */
template<class T>
void delete_matrix(T** mat, const short m) {
    for (short i = 0; i < m; i ++) {
        delete[] mat[i];
    }
    delete[] mat;
}


/**
 * @brief Хранит координаты точки
 * 
 */
struct Point {
    double x;
    double y;

    Point(const double& x, const double& y)
        : x(x), y(y) {}
};


/**
 * @brief Хранит точку и вектор
 * 
 */
struct Ray {
    double x;
    double y;
    double vx;
    double vy;

    Ray(const double& x, const double& y,
        const double& vx, const double& vy)
            : x(x), y(y), vx(vx), vy(vy) {}
};


/**
 * @brief Бокс для сетки.
 *  
 * Хранит индекс бокса и число пересечений линий с боксом.
 */
struct LightBox {
    short i, j, counter;

    LightBox()
        : i(0), j(0), counter(0) {}

    LightBox(const short& i, const short& j, const short& counter)
        : i(i), j(j), counter(counter) {}
    
    void set(const short& i, const short& j, const short& counter) {
        this->i = i;
        this->j = j;
        this->counter = counter;
    }
};


/**
 * @brief Сетка
 * 
 * Считает пересечения с линиями, ищет топ пересекаемых боксов
 */
struct Grid {
public:
    // Левый верхний угол сетки (координаты)
    const Point p_min;
    const Point p_max;
    // Количество боксов на сетке по разным осям
    const short ny;
    const short nx;
    // Размер бокса
    const double xs;
    const double ys;
    // Счетчики пересечений
    short** boxes;

    /**
     * @brief Construct a new Grid object
     * 
     * @param p_min Левый верхний угол сетки (координаты)
     * @param nx Количество боксов на сетке по x
     * @param ny Количество боксов на сетке по y
     * @param xs Размер бокса по x
     * @param ys Размер бокса по y
     */
    Grid(const Point& p_min, const short& nx, const short& ny,
         const double& xs, const double& ys)
            : p_min(p_min), p_max(p_min.x + nx*xs, p_min.y + ny*ys),
              nx(nx), ny(ny), xs(xs), ys(ys),
              boxes(create_matrix<short>(nx, ny, 0)) 
    {
        // Для поддержки разных размеров, реализуйте функцию подсчета расстояния,
        //     которая не использует факт, что xs = ys.
        // По умолчанию, используется fast_dist, которая предполагает, что xs = ys.
        // Так же вы можете просто растянуть одну из осей и использовать равные значения.
        if (xs != ys) {
            throw std::invalid_argument("XS != YS");
        }
    }

    ~Grid() {
        delete_matrix<short>(boxes, nx);
    }

    /**
     * @brief Считает пересечение линии с границой сетки
     * 
     * Инвариантен относительно направления линии.
     * 
     * @param ray Линия
     * @return Ray* Вернет nullptr, если пересечения нет,
     *         иначе вернет луч, который начинается на границе сетки
     *         и направлен во внутренность сетки. 
     */
    Ray* intersection(const Ray& ray) {
        double x, y, t;

        if (ray.vx != 0.0) {
            const double vx_inv = 1 / ray.vx;
            
            // x_min border
            t = (p_min.x - ray.x) * vx_inv;
            y = ray.y + t*ray.vy;
            if (y >= p_min.y && y <= p_max.y) {
                x = p_min.x;
                if (ray.vx >= 0) {
                    return new Ray(x, y, ray.vx, ray.vy); 
                } else {
                    return new Ray(x, y, -ray.vx, -ray.vy);
                }
            }

            // x_max border
            t = (p_max.x - ray.x) * vx_inv;
            y = ray.y + t*ray.vy;
            if (y >= p_min.y && y <= p_max.y) {
                x = p_max.x;
                if (ray.vx >= 0) {
                    return new Ray(x, y, -ray.vx, -ray.vy); 
                } else {
                    return new Ray(x, y, ray.vx, ray.vy); 
                }
            }

        }

        if (ray.vy != 0.0) {
            const double vy_inv = 1 / ray.vy;

            // y_min border
            t = (p_min.y - ray.y) * vy_inv;
            x = ray.x + t*ray.vx;
            if (x >= p_min.x && x <= p_max.x) {
                y = p_min.y;
                if (ray.vy >= 0) {
                    return new Ray(x, y, ray.vx, ray.vy);
                } else {
                    return new Ray(x, y, -ray.vx, -ray.vy);
                }
            }

            // y_max border
            t = (p_max.y - ray.y) * vy_inv;
            x = ray.x + t*ray.vx;
            if (x >= p_min.x && x <= p_max.x) {
                y = p_max.y;
                if (ray.vy >= 0) {
                    return new Ray(x, y, -ray.vx, -ray.vy);
                } else {
                    return new Ray(x, y, ray.vx, ray.vy);
                }
            }
        }

        return nullptr;
    }

    /**
     * @brief Считает пересечения ОДНОЙ линии с боксами
     * 
     * @param ray Точка пересечения с границей сетки
     *            + направление движения во внутрь сетки.
     *            Если не на границе или не смотрит во внутрь,
     *            нужно передать сюда результат вызова intersection(ray)
     */
    void voxel_traversal(const Ray& ray) {
        // Индексы первого пересеченного бокса
        short i = short((ray.x - p_min.x) / xs);
        short j = short((ray.y - p_min.y) / ys);

        // Проверяем, что не вышли за пределы сетки
        // Если вышли, корректируем
        i = to_safe(i, nx);
        j = to_safe(j, ny);

        // x variables init
        short stepX;
        double tDeltaX;
        double tMaxX;
        if (ray.vx > 0.0) {
            stepX = 1;
            tDeltaX = xs / ray.vx;
            tMaxX = (p_min.x + (i + 1) * xs - ray.x) / ray.vx;
        } else if (ray.vx < 0.0) {
            stepX = -1;
            tDeltaX = xs / -ray.vx;
            tMaxX = (p_min.x + i * xs - ray.x) / ray.vx;
        } else {
            stepX = 0;
            tDeltaX = MAX_TIME;
            tMaxX = MAX_TIME;
        }

        // y variables init
        short stepY;
        double tDeltaY;
        double tMaxY;
        if (ray.vy > 0.0) {
            stepY = 1;
            tDeltaY = ys / ray.vy;
            tMaxY = (p_min.y + (j + 1) * ys - ray.y) / ray.vy;
        } else if (ray.vy < 0.0) {
            stepY = -1;
            tDeltaY = ys / -ray.vy;
            tMaxY = (p_min.x + j * ys - ray.y) / ray.vy;
        } else {
            stepY = 0;
            tDeltaY = MAX_TIME;
            tMaxY = MAX_TIME;
        }

        // Обход вокселей(боксов)
        while (i >= 0 && i < nx && j >= 0 && j < ny) {
            boxes[i][j] += 1;
            if (tMaxX < tMaxY) {
                i += stepX;
                tMaxX += tDeltaX;
            } else {
                j += stepY;
                tMaxY += tDeltaY;
            }
        }
    }

    /**
     * @brief Считает пересечения МНОГИХ линий с боксами
     * 
     * @param rays Линии
     */
    void count_intersections(const std::vector<Ray> rays) {
        Ray* inter_ray;
        for (short i = 0; i < rays.size(); i++) {
            inter_ray = intersection(rays[i]);
            if (inter_ray != nullptr) {
                voxel_traversal(*inter_ray);    
            }
            delete inter_ray;
        }
    }

    /**
     * @brief Ищет топ k боксов по кол-ву пересечений с линиями
     * 
     * @param k 
     * @return std::vector<LightBox>* 
     */
    std::vector<LightBox>*
    top(const short& k) const {
        short i, j, r, l;
        short counter;
        std::vector<LightBox>& res = *(new std::vector<LightBox>);
        res.reserve(k);

        for (short r = 0; r < k; r++) {
            // HACK: -1 < boxes[i][j] for any i, j
            res.push_back(LightBox(-1000, -1000, -1));
        }
        
        const short k_minus_1 = k - 1;
        for (i = 0; i < nx; i++) {
            for(j = 0; j < ny; j++) {
                // add element to top_boxes
                counter = boxes[i][j];
                // делаем проверку вместо пробежки по k элементам
                if (counter >= res[k_minus_1].counter) {
                    for (l = 0; l < k; l++) {
                        // Если нужно поместить элемент на место l
                        if (counter >= res[l].counter) {
                            // Смещаем элементы
                            for (r = k_minus_1; r > l; r--) {
                                res[r] = res[r - 1];
                            }
                            // Добавляем на место l
                            res[l].set(i, j, counter);
                            break;
                        }
                    }
                }
            }
        }

        return &res;
    }

    /**
     * @brief Быстро вычисляет расстояние(НЕ ЕВКЛИДОВО, но очень похоже)
     * 
     * Правильно работает только, если xs == ys
     * 
     * @param b1 
     * @param b2 
     * @return double 
     */
    double fast_dist(const LightBox& b1, const LightBox& b2) const {
        return std::pow(b1.i - b2.i, 2) + std::pow(b1.j - b2.j, 2);
    }

    /**
     * @brief Считает топ боксов, так чтобы fast_dist между ними
     *        было больше D
     * 
     * @param k_start Кол-во топ боксов на первом этапе
     * @param k_end Кол-во боксов после отфильтровки по расстоянию
     * @return std::vector<LightBox>* 
     */
    std::vector<LightBox>*
    long_distance_top(const short& k_start, const short& k_end) const {
        bool flag;
        std::vector<LightBox>* full_top = top(k_start);
        std::vector<LightBox>& res = *(new std::vector<LightBox>);
        res.reserve(k_end);

        for(int i = 0; res.size() < k_end && i < k_start; i++) {
            flag = true;
            for(int j = 0; j < res.size(); j++) {
                if (fast_dist(res[j], (*full_top)[i]) <= D) {
                    flag = false;
                }
            }
            if (flag) {
                res.push_back((*full_top)[i]);
            }
        }

        delete full_top;
        return &res;
    }

private:
    short to_safe(short ind, short n) {
        if (ind >= n)
            return n - 1;
        if (ind < 0)
            return 0;
        return ind;
    }
};


/**
 * @brief Для тестирования промежуточных результатов
 * 
 * @param rays 
 */
void test_voxel_traversal(const std::vector<Ray>& rays) {
    Grid grid(Point(PMIN_X, PMIN_Y),
              NBOXES_X, NBOXES_Y,
              XS, YS);
    grid.count_intersections(rays);

    
    std::ofstream ofs = std::ofstream();
    ofs.open(TESTING_OUTPUT_NAME);
    for (short i = 0; i < NBOXES_X; i++)
        for (short j = 0; j < NBOXES_Y; j++)
            ofs << i << ' ' << j << ' ' << grid.boxes[i][j] << '\n';
    ofs.close();      
}


/**
 * @brief Алгоритм поиска точек
 * 
 * @param rays 
 * @param M 
 * @param N 
 */
void search_points(const std::vector<Ray>& rays, const short& M, const short& N) {
    Grid grid(Point(PMIN_X, PMIN_Y),
              NBOXES_X, NBOXES_Y,
              XS, YS);
    grid.count_intersections(rays);
    short k = short(N);
    std::vector<LightBox>* top = grid.long_distance_top(short(TOP_1*k), k);
    
    std::ofstream ofs = std::ofstream();
    ofs.open(OUTPUT_NAME);
    ofs << N << std::endl;
    for (short i = 0; i < k; i++) {
        LightBox box = (*top)[i];
        // ???: Нужно ли выводить counter в файл?
        ofs << box.i*XS + PMIN_X + XMIN << ' ' << box.j*YS + PMIN_Y + YMIN << ' ' << box.counter << std::endl;
    }
    ofs.close();
    delete top;
}


/**
 * @brief Загрузка линий из файла
 * 
 * @param filename 
 * @return std::tuple<short, short, std::vector<Ray>*> 
 */
std::tuple<short, short, std::vector<Ray>*>
load_lines_2(const std::string filename) {
    short M;  // кол-во прямых
    short N;  // кол-во точек
    std::vector<Ray>* a = new std::vector<Ray>();  // вектор из линий
    double x1, y1, x2, y2;

    std::ifstream ifs;
    ifs.open(filename);
    ifs >> M >> N;
    a->reserve(M);
    for(short i = 0; i < M; i++) {
        ifs >> x1 >> y1 >> x2 >> y2;
        // TODO: magick constants
        a->push_back(Ray(x1 - XMIN, y1 - YMIN, x2 - x1, y2 - y1));
    }
    ifs.close();
    
    return {M, N, a};
}


int main() {
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Start\n";
    auto [M, N, a] = load_lines_2(INPUT_NAME);
    // test_voxel_traversal(*a);
    search_points(*a, M, N);
    delete a;
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "End. Time of execution: "
        <<  std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << "ms\n";
}