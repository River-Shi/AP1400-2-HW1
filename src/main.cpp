#include <iostream>
#include <gtest/gtest.h>
#include "hw1.h"

int main(int argc, char **argv) {
    if (true) // make false to run unit-tests
    {
        Matrix a = algebra::random(3, 3, -10, 10);
        algebra::show(a);
    } else {
        ::testing::InitGoogleTest(&argc, argv);
        std::cout << "RUNNING TESTS ..." << std::endl;
        int ret{RUN_ALL_TESTS()};
        if (!ret)
            std::cout << "<<<SUCCESS>>>" << std::endl;
        else
            std::cout << "FAILED" << std::endl;
    }
    return 0;
}
