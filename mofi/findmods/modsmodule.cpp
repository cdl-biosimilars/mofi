#include "Python.h"
#include <vector>
#include "modifications.hpp"

const char* DOCSTRING = 
    "Finds all possible modifications for the given mass\n"
    "\n"
    "Arguments:\n"
    "mods -- A list of mod tuples with the format (mass, max amount)\n"
    "target_mass -- The mass to find modifications for\n"
    "massrange -- The tolerance of the mass value\n";

extern "C" {

    static PyObject* examine_modifications(PyObject* self, PyObject* args)
    {
        // Argument parsing
        PyObject* mods_list;
        double target_mass;
        double massrange = 5.0;
        if (!PyArg_ParseTuple(args, "O!d|d", &PyList_Type, &mods_list, &target_mass, &massrange)) {
            return NULL;
        }

        // Convert mods list to vector<Modification>
        size_t mods_list_size = PyList_Size(mods_list);
        std::vector<Modification> mods;
        for (size_t list_i = 0; list_i < mods_list_size; list_i++) {
            PyObject* mod_tuple = PyList_GetItem(mods_list, list_i);
            if (!PyTuple_Check(mod_tuple)) {
                PyErr_SetString(PyExc_TypeError, "Mod is not a tuple");
                return NULL;
            }
            PyObject* mass_item = PyTuple_GetItem(mod_tuple, 0);
            if (!mass_item) {
                PyErr_SetString(PyExc_IndexError, "Mod mass missing");
                return NULL;
            }
            PyObject* max_item = PyTuple_GetItem(mod_tuple, 1);
            if (!max_item) {
                PyErr_SetString(PyExc_IndexError, "Mod max missing");
                return NULL;
            }
            double mass = PyFloat_AsDouble(mass_item);
            if (PyErr_Occurred()) {
                return NULL;
            }
            long max = PyLong_AsLong(max_item);
            if (PyErr_Occurred()) {
                return NULL;
            }
            Modification mod = {mass, max};
            mods.push_back(mod);
        }

        // Find all the solutions
        std::vector<mod_state_t> solutions;
        solutions = find_modifications(target_mass, mods, massrange);

        // Convert vector back to list
        PyObject* res_list;
        res_list = PyList_New(solutions.size());
        for (size_t sol_i = 0; sol_i < solutions.size(); sol_i++) {
            const mod_state_t& sol = solutions[sol_i];
            PyObject* sol_list = PyList_New(sol.size());
            for (size_t count_i = 0; count_i < sol.size(); count_i++) {
                PyObject* count = PyLong_FromLong(sol[count_i]);
                PyList_SetItem(sol_list, count_i, count);
            }
            PyList_SetItem(res_list, sol_i, sol_list);
        }

        return res_list;
    }

    static PyMethodDef module_methods[] = {
        {"examine_modifications", (PyCFunction)examine_modifications, METH_VARARGS, DOCSTRING},
        {NULL, NULL, 0, NULL}
    };

    static struct PyModuleDef findmodsmodule = {
        PyModuleDef_HEAD_INIT,
        "findmods",
        "A C++ module that implements a fast combinatorial search",
        -1,
        module_methods
    };

    PyMODINIT_FUNC PyInit_findmods(void)
    {
        return PyModule_Create(&findmodsmodule);
    }

}
