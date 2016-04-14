/* Copyright (C) 2014, 2016 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <iostream>             // FIXME: this is just for error messages.
#include <vector>
#include <cassert>              // assert
#include <cstdlib>              // free

#include "ConfigJSON.hh"
#include "error_handling.hh"
#include "pism_utilities.hh"
#include "io/PIO.hh"

namespace pism {

/*! Given a 'path' "alice.bob.charlie", look for the JSON object 'bob'
 * containing a key 'charlie' in the object 'alice'. Use the 'object'
 * argument as the root.
 *
 * In other words, 'path' describes a node of a tree with 'object' as
 * the root, and this function returns the pointer to the node if
 * found, and NULL otherwise.
 */
static json_t* find_value(json_t *object,
                          const std::string &path) {
  if (object == NULL) {
    return NULL;
  }

  std::vector<std::string> split_path = split(path, '.');
  std::vector<std::string>::const_iterator j = split_path.begin();
  json_t *current_object = object;
  std::string current_path = *j;

  while (j != split_path.end()) {
    if (j != split_path.begin()) {
      current_path += "." + (*j);
    }

    current_object = json_object_get(current_object, j->c_str());
    if (current_object == NULL) {
      throw RuntimeError::formatted("ERROR: cannot find %s", current_path.c_str());
    }
    ++j;
  }

  return current_object;
}

ConfigJSON::ConfigJSON(units::System::Ptr unit_system)
  : Config(unit_system) {
  m_data = NULL;
}

ConfigJSON::~ConfigJSON() {
  json_decref(m_data);
}

/*! Initialize the database by reading from a file 'filename'.
 */
void ConfigJSON::init_from_file(const std::string &filename) {

  // free existing data if present
  if (m_data != NULL) {
    json_decref(m_data);
  }

  json_error_t error;
  m_data = json_load_file(filename.c_str(), JSON_DECODE_INT_AS_REAL, &error);

  if (m_data == NULL) {
    throw RuntimeError::formatted("Error loading config from '%s'"
                                  " at line %d, column %d.",
                                  filename.c_str(), error.line, error.column);
  }
}

/*! Initialize the database using a string 'string'.
 */
void ConfigJSON::init_from_string(const std::string &string) {
  // free existing data if present
  if (m_data != NULL) {
    json_decref(m_data);
  }

  json_error_t error;
  m_data = json_loads(string.c_str(), JSON_DECODE_INT_AS_REAL, &error);

  if (m_data == NULL) {
    throw RuntimeError::formatted("Error loading config from '%s'"
                                  " at line %d, column %d.",
                                  string.c_str(), error.line, error.column);
  }
}

/*! Return the JSON string representation of the configuration database.
 */
std::string ConfigJSON::dump() const {
  if (m_data == NULL) {
    return "";
  }

  char *tmp = json_dumps(m_data, JSON_INDENT(2) | JSON_ENSURE_ASCII | JSON_SORT_KEYS);
  std::string result;
  if (tmp != NULL) {
    result = tmp;
    free(tmp);
  }
  return result;
}

void ConfigJSON::read_impl(const PIO &nc) {
  std::string config_string = nc.get_att_text("PISM_GLOBAL", "pism_config");
  this->init_from_string(config_string);
}

void ConfigJSON::write_impl(const PIO &nc) const {
  std::string config_string = this->dump();

  nc.put_att_text("PISM_GLOBAL", "pism_config", config_string);
}

bool ConfigJSON::is_set_impl(const std::string &name) const {

}

Config::Doubles ConfigJSON::all_doubles_impl() const {

}

Config::Strings ConfigJSON::all_strings_impl() const {

}

Config::Booleans ConfigJSON::all_booleans_impl() const {

}

/*! Store a 'value' corresponding to the key 'name' in the database.
 *
 * If a name refers to an object "alice.bob", the object "alice" must
 * exist in the tree already, but "bob" may not exist before this call
 * and will be created. In other words, this method allows adding new
 * leaves only.
 */
void ConfigJSON::set_double_impl(const std::string &name, double value) {
  if (m_data == NULL) {
    return;
  }

  std::vector<std::string> path = split(name, '.');
  if (path.size() == 0) {
    // stop if 'name' is empty
    return;
  }

  std::string key = path.back();
  path.pop_back();

  json_t *object = NULL;
  if (path.empty() == true) {
    object = m_data;
  } else {
    object = find_value(m_data, join(path, "."));
  }

  if (object != NULL) {
    json_object_set_new(object, key.c_str(), json_pack("f", value));
  }
}

/*! Store a 'value' corresponding to the key 'name' in the database.
 *
 * If a name refers to an object "alice.bob", the object "alice" must
 * exist in the tree already, but "bob" may not exist before this call
 * and will be created. In other words, this method allows adding new
 * leaves only.
 */
void ConfigJSON::set_boolean_impl(const std::string &name, bool value) {
  if (m_data == NULL) {
    return;
  }

  std::vector<std::string> path = split(name, '.');
  if (path.size() == 0) {
    // stop if 'name' is empty
    return;
  }

  std::string key = path.back();
  path.pop_back();

  json_t *object = NULL;
  if (path.empty() == true) {
    object = m_data;
  } else {
    object = find_value(m_data, join(path, "."));
  }

  if (object != NULL) {
    json_object_set_new(object, key.c_str(), json_pack("b", value));
  }
}

/*! Store a 'value' corresponding to the key 'name' in the database.
 *
 * If a name refers to an object "alice.bob", the object "alice" must
 * exist in the tree already, but "bob" may not exist before this call
 * and will be created. In other words, this method allows adding new
 * leaves only.
 */
void ConfigJSON::set_string_impl(const std::string &name, const std::string &value) {
  if (m_data == NULL) {
    return;
  }

  std::vector<std::string> path = split(name, '.');
  if (path.size() == 0) {
    // stop if 'name' is empty
    return;
  }

  std::string key = path.back();
  path.pop_back();

  json_t *object = NULL;
  if (path.empty() == true) {
    object = m_data;
  } else {
    object = find_value(m_data, join(path, "."));
  }

  if (object != NULL) {
    json_object_set_new(object, key.c_str(), json_pack("s", value.c_str()));
  }
}

/*! Get the real number corresponding to the key 'name'.
 */
double ConfigJSON::get_double_impl(const std::string &name) const {
  json_t *value = find_value(m_data, name);
  if (value == NULL) {
    throw RuntimeError::formatted("%s was not found", name.c_str());
  }

  if (json_is_number(value) != 0) {
    double number = 0.0;
    if (json_unpack(value, "F", &number) == 0) {
      return number;
    } else {
      throw RuntimeError::formatted("failed to convert %s to double", name.c_str());
    }
  } else {
    throw RuntimeError::formatted("%s is not a number", name.c_str());
  }
}

/*! Get the string corresponding to the key 'name'.
 */
std::string ConfigJSON::get_string_impl(const std::string &name) const {
  json_t *value = find_value(m_data, name);
  if (value == NULL) {
    throw RuntimeError::formatted("%s was not found", name.c_str());
  }

  if (json_is_string(value) != 0) {
    const char *tmp_string = "";
    if (json_unpack(value, "s", &tmp_string) == 0) {
      return std::string(tmp_string);
    } else {
      throw RuntimeError::formatted("failed to convert %s to string", name.c_str());
    }
  } else {
    throw RuntimeError::formatted("%s is not a string", name.c_str());
  }
}

/*! Get the Boolean corresponding to the key 'name'.
 */
bool ConfigJSON::get_boolean_impl(const std::string &name) const {
  json_t *value = find_value(m_data, name);
  if (value == NULL) {
    throw RuntimeError::formatted("%s was not found", name.c_str());
  }

  if (json_is_boolean(value) != 0) {
    bool boolean = false;
    if (json_unpack(value, "b", &boolean) == 0) {
      return boolean;
    } else {
      throw RuntimeError::formatted("failed to convert %s to bool", name.c_str());
    }
  } else {
    throw RuntimeError::formatted("%s is not a boolean", name.c_str());
  }
}

} // end of namespace pism
