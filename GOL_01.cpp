#include <algorithm>
#include <array>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <thread>
#include <mutex>


class out_of_grid : public std::out_of_range
{
public:
    const int x, y, size_x, size_y;

    out_of_grid(int x, int y, int size_x, int size_y)
        :
        std::out_of_range(
            std::string("Grid index (") + std::to_string(x) + ", " + std::to_string(y) + ") is out of bounds for Grid of size[" + std::to_string(size_x) + ", " + std::to_string(size_y) + "]"
        ),
        x(x),
        y(y),
        size_x(size_x),
        size_y(size_y)
    {}
};


// Klasa Grid jest kontenerem komórek i planszą na której rozgrywa się "gra" (automat komórkowy). Oparta o wektor, logicznie implementuje 2-wymiarową tablicę o rozmiarze ustalanym w momencie konstrukcji.

template <typename CellT>
class Grid
{
public:
    const size_t size_x;
    const size_t size_y;

protected:
    std::vector<CellT> _storage;

    size_t _storage_index(int x, int y) const {
        return y * size_x + x;
    }

public:
    Grid(size_t size_x, size_t size_y)
        :
        size_x{ size_x },
        size_y{ size_y },
        _storage(size_x* size_y)
    {}

    CellT get_cell_unchecked(int x, int y) const {
        return _storage[_storage_index(x, y)];
    }

    CellT get_cell(int x, int y) const {
        check_bounds(x, y);
        return get_cell_unchecked(x, y);
    }

    CellT get_cell(int x, int y, const CellT& default_value) const {
        if (!is_in_bounds(x, y)) {
            return default_value;
        }

        return get_cell_unchecked(x, y);
    }

    void set_cell_unchecked(int x, int y, const CellT& value) {
        _storage[_storage_index(x, y)] = value;
    }

    void set_cell(int x, int y, const CellT& value) {
        check_bounds(x, y);
        set_cell_unchecked(x, y, value);
    }

    bool is_in_bounds(int x, int y) const {
        return 0 <= x && x < static_cast<int>(size_x) && 0 <= y && y < static_cast<int>(size_y);
    }

    void check_bounds(int x, int y) const {
        if (!is_in_bounds(x, y)) {
            throw out_of_grid(x, y, size_x, size_y);
        }
    }
};


// Klasy MooreNeighbourhood, VonNeumannNeighbourhood implementują różne rodzaje sąsiedztwa

class MooreNeighbourhood
{
public:
    constexpr static auto get_neighbourhood(int x, int y)
    {
        // Przy włączonym standardzie C++20 (std:c++latest w Visual Studio) konstruktor std::array przyjmujący listę inicjalizacyjną jest w stanie sam wydedukować liczbę elementów
        return std::array{
            std::make_tuple(x - 1, y - 1),
            std::make_tuple(x - 1, y),
            std::make_tuple(x - 1, y + 1),
            std::make_tuple(x,   y - 1),
            std::make_tuple(x,   y + 1),
            std::make_tuple(x + 1, y - 1),
            std::make_tuple(x + 1, y),
            std::make_tuple(x + 1, y + 1)
        };
    }
};


class VonNeumannNeighbourhood
{
public:
    constexpr static auto get_neighbourhood(int x, int y)
    {
        return std::array{
            std::make_tuple(x - 1, y),
            std::make_tuple(x,   y - 1),
            std::make_tuple(x + 1, y),
            std::make_tuple(x,   y + 1)
        };
    }
};


// Klasy WallBoundaryCondition, LoopBoundaryCondition zaimplementowane zostały jako pośrednicy, przez które kierowane są zapytania o wartosć sąsiadującej komórki podczas wykonywania kroku symulacji.

class WallBoundaryCondition
{
public:
    template<typename CellT>
    constexpr static auto get_cell(const Grid<CellT>& grid, int x, int y)
    {
        try {
            return grid.get_cell(x, y);
        }
        catch (const out_of_grid& exc) {
            return default_cell_value<CellT>();
        }
    }

    // Poprzez wydzielenie statycznej metody default_cell_value, możliwa jest zmiana tego aspektu warunku brzegowego "ściana" poprzez użycie polimorfizmu statycznego
    template<typename CellT>
    constexpr static auto default_cell_value() {
        return CellT();
    }
};


class LoopBoundaryCondition
{
    constexpr static int mod(int k, int n) {
        return ((k %= n) < 0) ? k + n : k;
    }

public:
    template<typename CellT>
    constexpr static auto get_cell(const Grid<CellT>& grid, int x, int y)
    {
        x = mod(x, grid.size_x);
        y = mod(y, grid.size_y);

        return grid.get_cell(x, y);
    }
};


// Klasa GameOfLifeRule implementuję logikę automatu komórkowego "Gra w życie". Jest parametryzowana typem komórki (w wypadku tego automatu sensownym wyborem jest typ bool),
// klasą implementującą sąsiedztwo oraz klasą implementującą warunki brzegowe
template <
    typename CellT,
    typename NeighbourhoodT,
    typename BoundaryConditionT
>
class GameOfLifeRule
{
public:
    const unsigned int survive_min, survive_max, birth_min, birth_max;

    GameOfLifeRule(
        unsigned int survive_min,
        unsigned int survive_max,
        unsigned int birth_min,
        unsigned int birth_max
    )
        :
        survive_min{ survive_min },
        survive_max{ survive_max },
        birth_min{ birth_min },
        birth_max{ birth_max }
    {}

    CellT operator() (const Grid<CellT>& grid, int x, int y) const
    {
        CellT current_cell = grid.get_cell(x, y);

        auto neighbourhood = NeighbourhoodT::get_neighbourhood(x, y);
        unsigned int alive_count = std::count_if(
            begin(neighbourhood),
            end(neighbourhood),
            // Usage of lambda function capturing 'grid' by reference
            [&grid](const auto& p) {
                auto [n_x, n_y] = p;
                return bool(BoundaryConditionT::get_cell(grid, n_x, n_y));
            });

        if (current_cell) {
            return CellT(survive_min <= alive_count && alive_count <= survive_max);
        }
        else {
            return CellT(birth_min <= alive_count && alive_count <= birth_max);
        }
    }
};


// CellularAutomaton2d clas implements a celullar automaton.
// To avoid allocating new Grid object in each simulation step, a double-buffering technique have been utilized
template <typename CellT, typename RuleT>
class CellularAutomaton2d
{
public:
    const size_t size_x;
    const size_t size_y;

protected:
    const RuleT _rule;
    std::unique_ptr<Grid<CellT>> _current_grid;
    std::unique_ptr<Grid<CellT>> _next_grid;

    void _flip_grids() {
        std::swap(_current_grid, _next_grid);
    }

public:
    CellularAutomaton2d(size_t size_x, size_t size_y, RuleT rule)
        :
        size_x{ size_x },
        size_y{ size_y },
        _rule(rule),
        _current_grid{ std::make_unique<Grid<CellT>>(size_x, size_y) },
        _next_grid{ std::make_unique<Grid<CellT>>(size_x, size_y) }
    {}

    void iterate(unsigned int steps) {
        while (steps--) {
            single_step();
        }
    }

    void single_step() {
        std::mutex mtx;
        std::vector<std::thread> threads;
        threads.reserve(size_x);

        for (size_t x = 0; x < size_x; ++x) {
            auto thread_function = [this, &mtx](int x) {
                std::vector<CellT> new_row;
                new_row.reserve(size_y);

                for (size_t y = 0; y < size_y; ++y) {
                    new_row.emplace_back(_rule(*_current_grid, x, y));
                }

                {
                    std::unique_lock<std::mutex> lock(mtx);
                    for (size_t y = 0; y < size_y; ++y)
                    {
                        _next_grid->set_cell(x, y, new_row[y]);
                    }
                }
            };
            threads.emplace_back(thread_function, x);
        }

        for (std::thread& thread : threads) {
            thread.join();
        }

        _flip_grids();
    }

    Grid<CellT>& get_current_grid() {
        return *_current_grid;
    }

    const Grid<CellT>& get_current_grid() const {
        return *_current_grid;
    }

    CellT get_cell(int x, int y) const {
        return _current_grid->get_cell(x, y);
    }

    CellT get_cell(int x, int y, const CellT& default_value) const {
        return _current_grid->get_cell(x, y, default_value);
    }

    void set_cell(int x, int y, const CellT& value) {
        _current_grid->set_cell(x, y, value);
    }

    template<typename NumberPairsIterable>
    void set_shape(const NumberPairsIterable& shape, const CellT& value) {
        // Usage of <algorithm> library
        std::for_each(
            begin(shape),
            end(shape),
            // Usage of lambda capturing 'this' (pointer to the current object) by value (since it's a pointer, capturing it by reference is not necessary) and 'value' function argument by reference
            [this, &value](const auto& p) {
                auto [x, y] = p;
                this->set_cell(x, y, value);
            });
    }
};


// Template function
template <typename GridT>
void display_grid(const GridT& grid)
{
    for (size_t y = 0; y < grid.size_y; ++y) {
        for (size_t x = 0; x < grid.size_x; ++x) {
            std::cout << (grid.get_cell(x, y) ? 'X' : '.');
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


int main()
{
    using GolRule = GameOfLifeRule<bool, MooreNeighbourhood, WallBoundaryCondition>;

    auto conway_gol_rule = GolRule(2, 3, 3, 3);
    CellularAutomaton2d<bool, GolRule> gol(40, 10, conway_gol_rule);

    // Initialization
    std::tuple<int, int> glider[] = {
        {1, 0},
        {2, 1},
        {2, 2},
        {1, 2},
        {0, 2}
    };
    gol.set_shape(glider, true);

    display_grid(gol.get_current_grid());

    // // Iteration
    // automaton.iterate(10);

    // Iteration and Display
    for (size_t i = 0; i < 30; ++i) {
        gol.single_step();
        display_grid(gol.get_current_grid());
    }

    return 0;
}
