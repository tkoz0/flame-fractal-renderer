#include "json.hpp"

namespace tkoz::flame
{

Json::Json()
{
}

Json::Json(const Json& json):
    nlohmann::json(static_cast<const nlohmann::json&>(json))
{
}

Json::Json(std::istream& input):
    nlohmann::json(nlohmann::json::parse(input,nullptr,true,true))
{
}

Json::Json(const std::string& input):
    nlohmann::json(nlohmann::json::parse(input,nullptr,true,true))
{
}

Json::Json(const nlohmann::json& input): nlohmann::json(input)
{
}

[[nodiscard]] bool Json::isNull() const noexcept
{
    return this->is_null();
}

[[nodiscard]] bool Json::isBool() const noexcept
{
    return this->is_boolean();
}

[[nodiscard]] bool Json::isInt() const noexcept
{
    return this->is_number_integer();
}

[[nodiscard]] bool Json::isFloat() const noexcept
{
    return this->is_number();
}

[[nodiscard]] bool Json::isString() const noexcept
{
    return this->is_string();
}

[[nodiscard]] bool Json::isArray() const noexcept
{
    return this->is_array();
}

[[nodiscard]] bool Json::isObject() const noexcept
{
    return this->is_object();
}

[[nodiscard]] size_t Json::size() const
{
    if (isArray() || isObject())
        return nlohmann::json::size();
    else
        throw JsonError("Json::size(): type is not an array or object");
}

[[nodiscard]] JsonBool Json::boolValue() const
{
    if (isBool())
        return this->get<JsonBool>();
    else
        throw JsonError("Json::boolValue(): type is not bool");
}

[[nodiscard]] JsonInt Json::intValue() const
{
    if (isInt())
        return this->get<JsonInt>();
    else
        throw JsonError("Json::intValue(): type is not int");
}

[[nodiscard]] JsonFloat Json::floatValue() const
{
    if (isFloat())
        return this->get<JsonFloat>();
    else
        throw JsonError("Json::floatValue(): type is not float");
}

[[nodiscard]] JsonString Json::stringValue() const
{
    if (isString())
        return this->get<JsonString>();
    else
        throw JsonError("Json::stringValue(): type is not string");
}

[[nodiscard]] JsonArray Json::arrayValue() const
{
    if (isArray())
    {
        auto tmp1 = this->get<std::vector<nlohmann::json>>();
        JsonArray tmp2;
        tmp2.reserve(tmp1.size());
        for (size_t i = 0; i < tmp1.size(); ++i)
            tmp2.push_back(Json(tmp1[i]));
        return tmp2;
    }
    else
        throw JsonError("Json::arrayValue(): type is not array");
}

[[nodiscard]] JsonObject Json::objectValue() const
{
    if (isObject())
    {
        auto tmp1 = this->get<std::unordered_map<std::string,nlohmann::json>>();
        JsonObject tmp2;
        tmp2.reserve(tmp1.size());
        for (auto iter = tmp1.begin(); iter != tmp1.end(); ++iter)
            tmp2.insert(std::make_pair(iter->first,Json(iter->second)));
        return tmp2;
    }
    else
        throw JsonError("Json::objectValue(): type is not object");
}

bool Json::valueAt(size_t index, Json& value) const noexcept
{
    if (!isArray() || index >= size())
        return false;
    value = static_cast<Json>(this->at(index));
    return true;
}

bool Json::valueAt(const std::string& key, Json& value) const noexcept
{
    if (!isObject() || this->find(key) == this->end())
        return false;
    value = static_cast<Json>(this->at(key));
    return true;
}

bool Json::valueAt(const char *key, Json& value) const noexcept
{
    if (!isObject() || this->find(key) == this->end())
        return false;
    value = static_cast<Json>(this->at(key));
    return true;
}

[[nodiscard]] Json Json::operator[](size_t index) const
{
    if (!isArray())
        throw JsonError("Json::operator[](size_t): not an array");
    else if (index >= size())
        throw JsonError("Json::operator[](size_t): array index out of range: "
            + std::to_string(index));
    else
        return static_cast<Json>(this->at(index));
}

[[nodiscard]] Json Json::operator[](const std::string& key) const
{
    if (!isObject())
        throw JsonError("Json::operator[](std::string): not an object");
    else if (this->find(key) == this->end())
        throw JsonError("Json::operator[](std::string): key does not exist: "
            + key);
    else
        return static_cast<Json>(this->at(key));
}

[[nodiscard]] Json Json::operator[](const char *key) const
{
    if (!isObject())
        throw JsonError("Json::operator[](char*): not an object");
    else if (this->find(key) == this->end())
        throw JsonError("Json::operator[](char*): key does not exist: "
            + std::string(key));
    else
        return static_cast<Json>(this->at(key));
}

Json& Json::operator=(const Json& a)
{
    static_cast<nlohmann::json&>(*this) = static_cast<const nlohmann::json&>(a);
    return *this;
}

bool Json::operator==(const Json& a)
{
    return static_cast<const nlohmann::json&>(*this)
        == static_cast<const nlohmann::json&>(a);
}

std::ostream& operator<<(std::ostream& os, const Json& json)
{
    os << static_cast<const nlohmann::json&>(json);
    return os;
}

} // namespace tkoz::flame
