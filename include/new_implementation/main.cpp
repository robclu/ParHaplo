#include "binary_container.hpp"
int main()
{
    haplo::BinaryContainer tiny(62);
    
    for (int i = 0; i < 10000000; ++i) {
        tiny.set(i % 41, 1);
    }
    std::cout << sizeof(tiny) << "\n";
    
    std::vector<haplo::detail::TinyContainer<haplo::byte>> tc;
    std::cout << sizeof(tc) << "\n";
}