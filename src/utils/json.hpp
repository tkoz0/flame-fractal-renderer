#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include "../../nlohmann/json.hpp"

namespace tkoz::flame
{

class JsonError : public std::exception
{
private:
    std::string msg;
public:
    JsonError() {}
    JsonError(const std::string& s): msg(s) {}
    JsonError(const char *s): msg(s) {}
    inline const char *what() const noexcept
    {
        return msg.c_str();
    }
};

class Json;

// JSON data types
// (null is not included and must be handled separately)
typedef std::vector<Json> JsonArray;
typedef std::unordered_map<std::string,Json> JsonObject;
typedef bool JsonBool;
typedef int64_t JsonInt;
typedef double JsonFloat;
typedef std::string JsonString;

/*
wrapper for nlohmann::json that can be compiled just once for this project,
because nlohmann's library adds a lot to the compile time
- supports comments in JSON
- read only / immutable
*/
class Json: private nlohmann::json
{
public:
    // create json object from input stream or string
    // create empty JSON
    Json();
    // copy constructor
    Json(const Json& json);
    // create JSON from input stream
    Json(std::istream& input);
    // create JSON from string
    Json(const std::string& input);
    // convert a nlohmann::json object into this type
    Json(const nlohmann::json& input);
    // is value null
    bool isNull() const noexcept;
    // is value true or false
    bool isBool() const noexcept;
    // is value an integer
    bool isInt() const noexcept;
    // is value a floating point
    bool isFloat() const noexcept;
    // is value a string
    bool isString() const noexcept;
    // is value an array
    bool isArray() const noexcept;
    // is value an object
    bool isObject() const noexcept;
    // returns the length of the array/object
    // throws an exception for other types
    size_t size() const;
    // returns boolean value, exception if not a boolean
    bool boolValue() const;
    // returns integer value, exception if not an integer
    JsonInt intValue() const;
    // returns floating point value, exception if not a floating point
    JsonFloat floatValue() const;
    // returns string value, exception if not a string
    std::string stringValue() const;
    // returns array value, exception if not an array
    JsonArray arrayValue() const;
    // returns object value, exception if not an object
    JsonObject objectValue() const;
    // if this is a long enough array, sets value and returns true
    bool valueAt(size_t index, Json& value) const noexcept;
    // if this is an object with the given key, sets value and returns true
    bool valueAt(const std::string& key, Json& value) const noexcept;
    // if this is an object with the given key, sets value and returns true
    bool valueAt(const char *key, Json& value) const noexcept;
    // access array index, exception if not an array or not long enough
    Json operator[](size_t index) const;
    // access object key, exception if not an object or does not have key
    Json operator[](const std::string& key) const;
    // access object key, exception if not an object or does not have key
    Json operator[](const char *key) const;
    // assignment operator
    Json& operator=(const Json& a);
    // compare equality of JSON objects
    bool operator==(const Json& a);
    // use nlohmann::json output stream operator
    friend std::ostream& operator<<(std::ostream& os, const Json& json);
};

} // namespace tkoz::flame
