#pragma once

#include <iostream>
#include <vector>
#include <unordered_map>

#include "../nlohmann/json.hpp"

class Json;

// type to use for JSON arrays
typedef std::vector<Json> JsonArray;

// type to use for JSON objects
typedef std::unordered_map<std::string,Json> JsonObject;

/*
interface for using nlohmann::json with the features needed for this project
supports comments in JSON
*/
class Json: private nlohmann::json
{
public:
    // create json object from input stream or string
    // create empty JSON
    Json();
    // create JSON from input stream
    Json(std::istream& input);
    // create JSON from string
    Json(const std::string& input);
    // convert a nlohmann::json object into this type
    Json(const nlohmann::json& input);
    // is value null
    bool isNull() const;
    // is value true or false
    bool isBool() const;
    // is value an integer
    bool isInt() const;
    // is value a floating point
    bool isFloat() const;
    // is value a string
    bool isString() const;
    // is value an array
    bool isArray() const;
    // is value an object
    bool isObject() const;
    // returns boolean value, exception if not a boolean
    bool boolValue() const;
    // returns integer value, exception if not an integer
    int64_t intValue() const;
    // returns floating point value, exception if not a floating point
    double floatValue() const;
    // returns string value, exception if not a string
    std::string stringValue() const;
    // returns array value, exception if not an array
    JsonArray arrayValue() const;
    // returns object value, exception if not an object
    JsonObject objectValue() const;
    // if this is a long enough array, sets value and returns true
    bool valueAt(size_t index, Json& value) const;
    // if this is an object with the given key, sets value and returns true
    bool valueAt(const std::string& key, Json& value) const;
    // if this is an object with the given key, sets value and returns true
    bool valueAt(const char *key, Json& value) const;
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
