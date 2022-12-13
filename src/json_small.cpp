#include "json_small.hpp"

Json::Json()
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

bool Json::isNull() const
{
    return this->is_null();
}

bool Json::isBool() const
{
    return this->is_boolean();
}

bool Json::isInt() const
{
    return this->is_number_integer();
}

bool Json::isFloat() const
{
    return this->is_number_float();
}

bool Json::isString() const
{
    return this->is_string();
}

bool Json::isArray() const
{
    return this->is_array();
}

bool Json::isObject() const
{
    return this->is_object();
}

bool Json::boolValue() const
{
    return this->get<bool>();
}

int64_t Json::intValue() const
{
    return this->get<int64_t>();
}

double Json::floatValue() const
{
    return this->get<double>();
}

std::string Json::stringValue() const
{
    return this->get<std::string>();
}

JsonArray Json::arrayValue() const
{
    auto tmp1 = this->get<std::vector<nlohmann::json>>();
    JsonArray tmp2;
    tmp2.reserve(tmp1.size());
    for (size_t i = 0; i < tmp1.size(); ++i)
        tmp2.push_back(Json(tmp1[i]));
    return tmp2;
}

JsonObject Json::objectValue() const
{
    auto tmp1 = this->get<std::unordered_map<std::string,nlohmann::json>>();
    JsonObject tmp2;
    tmp2.reserve(tmp1.size());
    for (auto iter = tmp1.begin(); iter != tmp1.end(); ++iter)
        tmp2.insert(std::make_pair(iter->first,Json(iter->second)));
    return tmp2;
}

bool Json::valueAt(size_t index, Json& value) const
{
    if (!this->is_array() || index >= this->size())
        return false;
    value = Json(this->at(index));
    return true;
}

bool Json::valueAt(const std::string& key, Json& value) const
{
    if (!this->is_object() || this->find(key) == this->end())
        return false;
    value = Json(this->at(key));
    return true;
}

bool Json::valueAt(const char *key, Json& value) const
{
    if (!this->is_object() || this->find(key) == this->end())
        return false;
    value = Json(this->at(key));
    return true;
}

Json Json::operator[](size_t index) const
{
    return Json(this->at(index));
}

Json Json::operator[](const std::string& key) const
{
    return Json(this->at(key));
}

Json Json::operator[](const char *key) const
{
    return Json(this->at(key));
}

Json& Json::operator=(const Json& a)
{
    static_cast<nlohmann::json&>(*this) = static_cast<const nlohmann::json&>(a);
    return *this;
}

bool Json::operator==(const Json& a)
{
    return static_cast<const nlohmann::json&>(*this) \
        == static_cast<const nlohmann::json&>(a);
}

std::ostream& operator<<(std::ostream& os, const Json& json)
{
    os << static_cast<const nlohmann::json&>(json);
    return os;
}
