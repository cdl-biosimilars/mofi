#include <vector>
#include <string>

typedef std::vector<int> mod_state_t;

/**
 * Struct to store modification details
 */
struct Modification
{
    /**
     * Mass of the modification
     */
    double mass;

    /**
     * Maximum amount of the modification on a protein
     */
    long max;
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
std::vector<mod_state_t> find_modifications(
    const double target_mass,
    std::vector<Modification> mods,
    const double massrange=5.0);
