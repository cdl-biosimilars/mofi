#include <vector>
#include <cmath>
#include <iostream>
#include <cassert>
#include <algorithm>

#include "modifications.hpp"


/**
 * Tolerance to prevent rounding errors from affecting the result
 */
const double TOLERANCE = 0.0000001;


/**
 * Recursive function to check all possible combinations of
 * modifications if they match a specific mass. Shouldn't be used
 * directly, use find_modifications instead.
 *
 * @param target_mass mass to find solutions for
 * @param mods modifications to use
 * @param solutions container to insert the found solutions into
 * @param used Already used modifications in lower recursion levels
 * @param massrange tolerance of the target mass
 * @param index mod to start at
 */
void next_mod(
    const double target_mass,
    const std::vector<Modification>& mods,
    std::vector<mod_state_t>& solutions,
    mod_state_t& used,
    const double massrange,
    const size_t index=0,
    double usedmass=0)
{
    long mod_max = mods[index].max;
    size_t next_index = index + 1;
    for (long x = 0; x <= mod_max; x++) {
        used[index] = x;
        double remaining = target_mass - usedmass;
        // Check if the mass matches
        if (std::abs(remaining) <= massrange) {
            solutions.push_back(used);
        } else if (next_index < mods.size() && remaining >= massrange) {
            next_mod(
                target_mass, mods, solutions, used, massrange,
                next_index, usedmass);
        }
        // Add current modification to total mass
        usedmass += mods[index].mass;
    }
    used[index] = 0;
}

std::vector<mod_state_t> find_modifications(
    const double target_mass,
    std::vector<Modification> mods,
    double massrange)
{
    std::vector<mod_state_t> solutions;
    mod_state_t used(mods.size(), 0);
    massrange += TOLERANCE;
    next_mod(target_mass, mods, solutions, used, massrange);
    return solutions;
}
