#include <vector>
#include <string>

typedef std::vector<int> mod_state_t;

// Struct to store modification details
struct Modification
{
    double mass;  // Mass of the modification
    long max;  // Maximum amount of the modification on a protein
};

// Struct to store the results of a modifications search
struct SearchResult
{
    unsigned long long int search_space_size;
    std::vector<mod_state_t> combs;
};

/**
 * Finds all possible combinations of modifications to match a specific
 * mass. Wrapper to create the initial state for next_mod.
 *
 * @param target_mass mass to find solutions for
 * @param mods modifications to use
 * @param massrange tolerance of the target mass
 * @return container with all found solutions
 */
SearchResult find_modifications(
    const double target_mass,
    std::vector<Modification> mods,
    const double massrange=5.0);