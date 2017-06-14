#include <vector>

#include "modifications.hpp"

int main()
{
    std::vector<Modification> mods = {
        {134.3, 6},
        {242.6, 13},
        {229.2, 10},
        {181.8, 16},
        {171.5, 8},
        {251.1, 12},
        {146.9, 11},
        {250.4, 8},
    };
    
    std::vector<mod_state_t> solutions;
    solutions = find_modifications(4299.2, mods);
    std::cout << solutions.size() << std::endl;
    return 0;
}
