#include "scippp/version.hpp"

namespace scippp {

static constexpr std::string_view VERSION = "local";

std::string_view getVersion()
{
    return VERSION;
}

}
