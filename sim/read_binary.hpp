// #ifndef READ_BINARY_HPP
// #define READ_BINARY_HPP

// template <typename T>
// requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
// char * as_writable_buffer(T & value) {
//     // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
//     return reinterpret_cast<char *>(&value);
// }

// template <typename T>
// requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
// char const * as_buffer(T const & value) {
//     // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
//     return reinterpret_cast<char const *>(&value);
// }

// template <typename T>
// requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
// T read_binary_value(std::istream & is) {
//     T value{};
//     is.read(as_writable_buffer(value), sizeof(value));
//     return value;
// }

// template <typename T>
// requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
// void write_binary_value(T value, std::ostream & os) {
//     //if (!(std::isnan(value))) {std::cout << value << '\n';}
//     os.write(as_buffer(value), sizeof(value));
// }

// #endif //READ_BINARY_HPP
