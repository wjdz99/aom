#include "examples/multilayer_metadata.h"

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "aom/aom_integer.h"
#include "examples/multilayer_metadata.h"

namespace libaom_examples {

namespace {

constexpr int kMaxNumSpatialLayers = 4;

// Removes comments and trailing spaces from the line.
void cleanup_line(std::string &line) {
  // Remove everything after the first '#'.
  std::size_t comment_pos = line.find('#');
  if (comment_pos != std::string::npos) {
    line.resize(comment_pos);
  }
  // Remove spaces at the end of the line.
  while (!line.empty() && line.back() == ' ') {
    line.resize(line.length() - 1);
  }
}

// Finds the indentation level of the line, and sets 'has_list_prefix' to true
// if the line has a '-' indicating a new item in a list.
void get_indent(const std::string &line, int *indent, bool *has_list_prefix) {
  *indent = 0;
  *has_list_prefix = 0;
  while (
      *indent < (int)line.length() &&
      (line[*indent] == ' ' || line[*indent] == '\t' || line[*indent] == '-')) {
    if (line[*indent] == '-') {
      *has_list_prefix = true;
    }
    ++(*indent);
  }
}

struct ParsedValue {
  enum class Type { kNone, kInt, kFloatingPoint };
  Type type = Type::kNone;

  long int_value = 0;
  double double_value = 0.0f;

  double ValueAsFloatingPoint(int line_idx) {
    if (type == Type::kNone) {
      fprintf(
          stderr,
          "No value found where floating point value was expected at line %d\n",
          line_idx);
      exit(EXIT_FAILURE);
    }
    return (type == Type::kFloatingPoint) ? double_value : (double)int_value;
  }

  long IntValueInRange(long min, long max, int line_idx) {
    switch (type) {
      case Type::kInt:
        if (int_value < min || int_value > max) {
          fprintf(stderr,
                  "Integer value %ld out of range [%ld, %ld] at line %d\n",
                  int_value, min, max, line_idx);
          exit(EXIT_FAILURE);
        }
        return int_value;
      case Type::kFloatingPoint:
        fprintf(stderr,
                "Floating point value found where integer was expected at line "
                "%d\n",
                line_idx);
        exit(EXIT_FAILURE);
      case Type::kNone:
      default:
        fprintf(stderr,
                "No value found where integer was expected at line %d\n",
                line_idx);
        exit(EXIT_FAILURE);
    }
  }
};

/*
 * Parses the next line from the file, skipping empty lines.
 * Returns false if the end of the file was reached, or if the line was indented
 * less than 'min_indent', meaning that parsing should go back to the previous
 * function in the stack.
 *
 * 'min_indent' is the minimum indentation expected for the next line.
 * 'is_list' must be true if the line is allowed to contain list items ('-').
 * 'indent' MUST be initialized to -1 before the first call, and is then set to
 * the indentation of the line.
 * 'has_list_prefix' is set to true if the line starts a new list item with '-'.
 * 'line_idx' is set to the index of the last line read.
 * 'field_name' is set to the field name if the line contains a colon, or to an
 * empty string otherwise.
 * 'value' is set to the value on the line if present.
 */
bool parse_line(std::fstream &file, int min_indent, bool is_list, int *indent,
                bool *has_list_prefix, int *line_idx, std::string *field_name,
                ParsedValue *value) {
  *field_name = "";
  *value = {};
  std::string line;
  std::fstream::pos_type prev_file_position;
  const int prev_indent = *indent;
  while (prev_file_position = file.tellg(), std::getline(file, line)) {
    cleanup_line(line);
    get_indent(line, indent, has_list_prefix);
    line = line.substr(*indent);  // skip indentation
    // If the line is indented less than 'min_indent', it belongs to the outer
    // object, and parsing should go back to the previous function in the stack.
    if (!line.empty() && *indent < min_indent) {
      // Undo reading the last line.
      if (!file.seekp(prev_file_position, std::ios::beg)) {
        fprintf(stderr, "Failed to seek to previous file position\n");
        exit(EXIT_FAILURE);
      }
      return false;
    }

    ++(*line_idx);
    if (line.empty()) continue;

    if (prev_indent >= 0 && prev_indent != *indent) {
      fprintf(stderr, "Error: Bad indentation at line %d\n", *line_idx);
      exit(EXIT_FAILURE);
    }
    if (*has_list_prefix && !is_list) {
      fprintf(stderr, "Error: Unexpected list item at line %d\n", *line_idx);
      exit(EXIT_FAILURE);
    }

    std::string value_str = line;
    size_t colon_pos = line.find(':');
    if (colon_pos != std::string::npos) {
      *field_name = line.substr(0, colon_pos);
      value_str = line.substr(colon_pos + 1);
    }
    if (!value_str.empty()) {
      char *endptr;
      if (line.find('.') != std::string::npos) {
        value->double_value = strtod(&line[colon_pos + 1], &endptr);
        value->type = ParsedValue::Type::kFloatingPoint;
        if (*endptr != '\0') {
          fprintf(stderr,
                  "Error: Failed to parse floating point value from '%s' at "
                  "line %d\n",
                  value_str.c_str(), *line_idx);
          exit(EXIT_FAILURE);
        }
      } else {
        value->int_value = strtol(&line[colon_pos + 1], &endptr, 10);
        value->type = ParsedValue::Type::kInt;
        if (*endptr != '\0') {
          fprintf(stderr,
                  "Error: Failed to parse integer from '%s' at line %d\n",
                  value_str.c_str(), *line_idx);
          exit(EXIT_FAILURE);
        }
      }
    }
    return true;
  }
  return false;  // Reached the end of the file.
}

template <typename T>
std::vector<T> parse_integer_list(std::fstream &file, int min_indent,
                                  int *line_idx) {
  bool has_list_prefix;
  int indent = -1;
  std::string field_name;
  ParsedValue value;
  std::vector<T> result;
  while (parse_line(file, min_indent, /*is_list=*/true, &indent,
                    &has_list_prefix, line_idx, &field_name, &value)) {
    if (!field_name.empty()) {
      fprintf(
          stderr,
          "Error: Unexpected field name '%s' at line %d, expected a number\n",
          field_name.c_str(), *line_idx);
      exit(EXIT_FAILURE);
    } else if (!has_list_prefix) {
      fprintf(stderr, "Error: Missing list prefix '-' at line %d\n", *line_idx);
      exit(EXIT_FAILURE);
    } else {
      result.push_back((T)value.IntValueInRange(
          (long)std::numeric_limits<T>::min(),
          (long)std::numeric_limits<T>::max(), *line_idx));
    }
  }
  return result;
}

template <typename T>
std::pair<T, bool> value_present(const T &v) {
  return std::make_pair(v, true);
}

ColorProperties parse_color_properties(std::fstream &file, int min_indent,
                                       int *line_idx) {
  bool has_list_prefix;
  int indent = -1;
  std::string field_name;
  ParsedValue value;
  ColorProperties color = {};
  while (parse_line(file, min_indent, /*is_list=*/false, &indent,
                    &has_list_prefix, line_idx, &field_name, &value)) {
    if (field_name == "color_range") {
      color.color_range =
          value.IntValueInRange(/*min=*/0, /*max=*/1, *line_idx);
    } else if (field_name == "color_primaries") {
      color.color_primaries =
          value.IntValueInRange(/*min=*/0, /*max=*/255, *line_idx);
    } else if (field_name == "transfer_characteristics") {
      color.transfer_characteristics =
          value.IntValueInRange(/*min=*/0, /*max=*/255, *line_idx);
    } else if (field_name == "matrix_coefficients") {
      color.matrix_coefficients =
          value.IntValueInRange(/*min=*/0, /*max=*/255, *line_idx);
    } else {
      fprintf(stderr, "Error: Unknown field '%s' at line %d\n",
              field_name.c_str(), *line_idx);
      exit(EXIT_FAILURE);
    }
  }
  return color;
}

AlphaInformation parse_multilayer_layer_alpha(std::fstream &file,
                                              int min_indent, int *line_idx) {
  bool has_list_prefix;
  int indent = -1;
  std::string field_name;
  ParsedValue value;
  AlphaInformation alpha_info = {};
  while (parse_line(file, min_indent, /*is_list=*/false, &indent,
                    &has_list_prefix, line_idx, &field_name, &value)) {
    if (field_name == "alpha_use_idc") {
      alpha_info.alpha_use_idc =
          (AlphaUse)value.IntValueInRange(/*min=*/0, /*max=*/7, *line_idx);
    } else if (field_name == "alpha_bit_depth") {
      alpha_info.alpha_bit_depth =
          value.IntValueInRange(/*min=*/8, /*max=*/15, *line_idx);
    } else if (field_name == "alpha_clip_idc") {
      alpha_info.alpha_clip_idc =
          value.IntValueInRange(/*min=*/0, /*max=*/3, *line_idx);
    } else if (field_name == "alpha_incr_flag") {
      alpha_info.alpha_incr_flag =
          value.IntValueInRange(/*min=*/0, /*max=*/1, *line_idx);
    } else if (field_name == "alpha_transparent_value") {
      // At this point we may not have parsed 'alpha_bit_depth' yet, so the
      // exact range is checked later.
      alpha_info.alpha_transparent_value = value.IntValueInRange(
          std::numeric_limits<uint16_t>::min(),
          std::numeric_limits<uint16_t>::max(), *line_idx);
    } else if (field_name == "alpha_opaque_value") {
      // At this point we may not have parsed 'alpha_bit_depth' yet, so the
      // exact range is checked later.
      alpha_info.alpha_opaque_value = value.IntValueInRange(
          std::numeric_limits<uint16_t>::min(),
          std::numeric_limits<uint16_t>::max(), *line_idx);
    } else if (field_name == "alpha_color_description") {
      alpha_info.alpha_color_description =
          value_present(parse_color_properties(file, indent, line_idx));
    } else if (field_name == "label_type_id") {
      alpha_info.label_type_id = parse_integer_list<uint16_t>(
          file, /*min_indent=*/indent + 1, line_idx);
    } else {
      fprintf(stderr, "Error: Unknown field '%s' at line %d\n",
              field_name.c_str(), *line_idx);
      exit(EXIT_FAILURE);
    }
  }

  // Validation.
  if (alpha_info.alpha_bit_depth == 0) {
    fprintf(stderr,
            "Error: alpha_bit_depth must be specified (in range [8, 15]) for "
            "alpha info\n");
    exit(EXIT_FAILURE);
  }
  const int alpha_max = (1 << alpha_info.alpha_bit_depth) - 1;
  if (alpha_info.alpha_transparent_value > alpha_max) {
    fprintf(stderr, "Error: alpha_transparent_value %d out of range [0, %d]\n",
            alpha_info.alpha_transparent_value, alpha_max);
    exit(EXIT_FAILURE);
  }
  if (alpha_info.alpha_opaque_value > alpha_max) {
    fprintf(stderr, "Error: alpha_opaque_value %d out of range [0, %d]\n",
            alpha_info.alpha_opaque_value, alpha_max);
    exit(EXIT_FAILURE);
  }
  if ((!alpha_info.label_type_id.empty()) &&
      (alpha_info.alpha_use_idc != ALPHA_SEGMENTATION)) {
    fprintf(stderr,
            "Error: label_type_id can only be set if alpha_use_idc is %d\n",
            ALPHA_SEGMENTATION);
    exit(EXIT_FAILURE);
  }
  const int alpha_range = (std::abs((int)alpha_info.alpha_opaque_value -
                                    alpha_info.alpha_transparent_value) +
                           1);
  if (!alpha_info.label_type_id.empty() &&
      (int)alpha_info.label_type_id.size() != alpha_range) {
    fprintf(stderr,
            "Error: if present, label_type_id size must be "
            "equal to the range of alpha values between "
            "alpha_transparent_value and alpha_opaque_value (expected "
            "%d values, found %d values)\n",
            alpha_range, (int)alpha_info.label_type_id.size());
    exit(EXIT_FAILURE);
  }
  if (alpha_info.alpha_color_description.second &&
      (alpha_info.alpha_use_idc != ALPHA_STRAIGHT)) {
    fprintf(stderr,
            "Error: alpha_color_description can only be set if alpha_use_idc "
            "is %d\n",
            ALPHA_STRAIGHT);
    exit(EXIT_FAILURE);
  }
  return alpha_info;
}

double depth_representation_element_to_double(
    const DepthRepresentationElement &e) {
  // Let x be a variable that is computed using four variables s, e, m, and n,
  // as follows: If e is greater than 0 and less than 127, x is set equal to
  // (−1)^s*2^(e−31) * (1+m÷2^n).
  // Otherwise (e is equal to 0), x is set equal to (−1)^s*2^−(30+n)*m.
  if (e.exponent > 0) {
    return (e.sign_flag ? -1 : 1) * std::pow(2.0, (int)e.exponent - 31) *
           (1 + (double)e.mantissa / ((int64_t)1 << e.mantissa_len));
  } else {
    return (e.sign_flag ? -1 : 1) * e.mantissa *
           std::pow(2.0, -30 + (int)e.mantissa_len);
  }
}

DepthRepresentationElement double_to_depth_representation_element(double f) {
  const double orig = f;
  if (f == 0.0) {
    return { 0, 0, 0, 1 };
  }
  const bool sign = f < 0.0;
  if (sign) {
    f = -f;
  }
  int exp = 0;
  if (f >= 1.0) {
    while (f >= 2.0) {
      ++exp;
      f /= 2;
    }
  } else {
    while (f < 1.0) {
      ++exp;
      f *= 2.0;
    }
    exp = -exp;
  }
  if ((exp + 31) <= 0 || (exp + 31) > 126) {
    fprintf(stderr,
            "Error: Floating point value %f too out of range (too large or too "
            "small)",
            orig);
    exit(EXIT_FAILURE);
  }
  assert(f >= 1.0 && f < 2.0);
  f -= 1.0;
  uint32_t mantissa = 0;
  uint16_t mantissa_len = 0;
  // Technically the specification allows up to 32 but
  // aom_wb_write_literal only supports writing up to 31 bits.
  constexpr uint16_t kMaxMantissaLen = 31;
  do {
    const int bit = (f >= 0.5);
    mantissa = (mantissa << 1) + bit;
    f -= bit * 0.5;
    ++mantissa_len;
    f *= 2.0;
  } while (mantissa_len < kMaxMantissaLen && f > 0.0);
  const DepthRepresentationElement element = { sign, (uint8_t)(exp + 31),
                                               mantissa, mantissa_len };
  return element;
}

DepthInformation parse_multilayer_layer_depth(std::fstream &file,
                                              int min_indent, int *line_idx) {
  bool has_list_prefix;
  int indent = -1;
  std::string field_name;
  ParsedValue value;
  DepthInformation depth_info = {};
  while (parse_line(file, min_indent, /*is_list=*/false, &indent,
                    &has_list_prefix, line_idx, &field_name, &value)) {
    if (field_name == "z_near") {
      depth_info.z_near = value_present(double_to_depth_representation_element(
          value.ValueAsFloatingPoint(*line_idx)));
    } else if (field_name == "z_far") {
      depth_info.z_far = value_present(double_to_depth_representation_element(
          value.ValueAsFloatingPoint(*line_idx)));
    } else if (field_name == "d_min") {
      depth_info.d_min = value_present(double_to_depth_representation_element(
          value.ValueAsFloatingPoint(*line_idx)));
    } else if (field_name == "d_max") {
      depth_info.d_max = value_present(double_to_depth_representation_element(
          value.ValueAsFloatingPoint(*line_idx)));
    } else if (field_name == "depth_representation_type") {
      depth_info.depth_representation_type =
          value.IntValueInRange(/*min=*/0, /*max=*/15, *line_idx);
    } else if (field_name == "disparity_ref_view_id") {
      depth_info.disparity_ref_view_id =
          value.IntValueInRange(/*min=*/0, /*max=*/3, *line_idx);
    } else if (field_name == "depth_nonlinear_precision") {
      depth_info.depth_nonlinear_precision =
          value.IntValueInRange(/*min=*/8, /*max=*/23, *line_idx);
    } else if (field_name == "depth_nonlinear_representation_model") {
      depth_info.depth_nonlinear_representation_model =
          parse_integer_list<uint32_t>(file,
                                       /*min_indent=*/indent + 1, line_idx);
    } else {
      fprintf(stderr, "Error: Unknown field '%s' at line %d\n",
              field_name.c_str(), *line_idx);
      exit(EXIT_FAILURE);
    }
  }

  // Validation.
  if (depth_info.depth_representation_type == 3 &&
      depth_info.depth_nonlinear_precision == 0) {
    fprintf(stderr,
            "Error: depth_nonlinear_precision must be specified (in range [8, "
            "23]) when "
            "depth_representation_type is 3\n");
    exit(EXIT_FAILURE);
  }
  if ((depth_info.depth_representation_type == 3) !=
      (!depth_info.depth_nonlinear_representation_model.empty())) {
    fprintf(stderr,
            "Error: depth_nonlinear_representation_model must be set if and "
            "only if depth_representation_type is 3\n");
    exit(EXIT_FAILURE);
  }
  const uint32_t depth_max = (1 << depth_info.depth_nonlinear_precision) - 1;
  for (uint32_t v : depth_info.depth_nonlinear_representation_model) {
    if (v > depth_max) {
      fprintf(stderr,
              "Error: depth_nonlinear_representation_model value %d out of "
              "range [0, %d]\n",
              v, depth_max);
      exit(EXIT_FAILURE);
    }
  }

  return depth_info;
}

std::vector<LayerMetadata> parse_multilayer_layer_metadata(std::fstream &file,
                                                           int min_indent,
                                                           int *line_idx) {
  bool has_list_prefix;
  int indent = -1;
  std::string field_name;
  ParsedValue value;
  std::vector<LayerMetadata> layers;
  bool layer_has_alpha = false;
  bool layer_has_depth = false;
  while (parse_line(file, min_indent, /*is_list=*/true, &indent,
                    &has_list_prefix, line_idx, &field_name, &value)) {
    if (has_list_prefix) {
      // Start of a new layer.
      if (layers.size() >= kMaxNumSpatialLayers) {
        fprintf(stderr,
                "Error: Too many layers at line %d, the maximum is %d\n",
                *line_idx, kMaxNumSpatialLayers);
        exit(EXIT_FAILURE);
      }

      // Validate the previous layer.
      if (!layers.empty()) {
        if (layer_has_alpha !=
            (layers.back().layer_type == MULTILAYER_LAYER_TYPE_ALPHA &&
             layers.back().layer_metadata_scope >= SCOPE_GLOBAL)) {
          fprintf(stderr,
                  "Error: alpha info must be set if and only if layer_type is "
                  "%d and layer_metadata_scpoe is >= %d\n",
                  MULTILAYER_LAYER_TYPE_ALPHA, SCOPE_GLOBAL);
          exit(EXIT_FAILURE);
        }
        if (layer_has_depth !=
            (layers.back().layer_type == MULTILAYER_LAYER_TYPE_DEPTH &&
             layers.back().layer_metadata_scope >= SCOPE_GLOBAL)) {
          fprintf(stderr,
                  "Error: depth info must be set if and only if layer_type is "
                  "%d and layer_metadata_scpoe is >= %d\n",
                  MULTILAYER_LAYER_TYPE_DEPTH, SCOPE_GLOBAL);
          exit(EXIT_FAILURE);
        }
      }
      if (layers.size() == 1 && layers.back().layer_color_description.second) {
        fprintf(stderr,
                "Error: layer_color_description cannot be specified for the "
                "first layer\n");
        exit(EXIT_FAILURE);
      }

      layers.push_back({});
      layer_has_alpha = false;
      layer_has_depth = false;
    }
    if (layers.empty()) {
      fprintf(stderr, "Error: Missing list prefix '-' at line %d\n", *line_idx);
      exit(EXIT_FAILURE);
    }

    LayerMetadata *layer = &layers.back();
    // Check if string starts with field name.
    if ((field_name == "layer_type")) {
      layer->layer_type =
          (LayerType)value.IntValueInRange(/*min=*/0, /*max=*/31, *line_idx);
    } else if ((field_name == "luma_plane_only_flag")) {
      layer->luma_plane_only_flag =
          value.IntValueInRange(/*min=*/0, /*max=*/1, *line_idx);
    } else if ((field_name == "layer_view_type")) {
      layer->layer_view_type = (MultilayerViewType)value.IntValueInRange(
          /*min=*/0, /*max=*/7, *line_idx);
    } else if ((field_name == "group_id")) {
      layer->group_id = value.IntValueInRange(/*min=*/0, /*max=*/3, *line_idx);
    } else if ((field_name == "layer_dependency_idc")) {
      layer->layer_dependency_idc =
          value.IntValueInRange(/*min=*/0, /*max=*/7, *line_idx);
    } else if ((field_name == "layer_metadata_scope")) {
      layer->layer_metadata_scope =
          (MultilayerMetadataScope)value.IntValueInRange(/*min=*/0, /*max=*/3,
                                                         *line_idx);
    } else if ((field_name == "layer_color_description")) {
      layer->layer_color_description =
          value_present(parse_color_properties(file, indent, line_idx));
    } else if ((field_name == "alpha")) {
      layer_has_alpha = true;
      layer->global_alpha_info =
          parse_multilayer_layer_alpha(file,
                                       /*min_indent=*/indent + 1, line_idx);
    } else if (field_name == "depth") {
      layer_has_depth = true;
      layer->global_depth_info =
          parse_multilayer_layer_depth(file,
                                       /*min_indent=*/indent + 1, line_idx);
      if ((layer->global_depth_info.d_min.second ||
           layer->global_depth_info.d_max.second) &&
          layer->global_depth_info.disparity_ref_view_id ==
              (layers.size() - 1)) {
        fprintf(stderr,
                "disparity_ref_view_id must be different from the layer's id "
                "for layer %d (zero-based index)\n",
                (int)layers.size() - 1);
        exit(EXIT_FAILURE);
      }
    } else {
      fprintf(stderr, "Error: Unknown field %s at line %d\n",
              field_name.c_str(), *line_idx);
      exit(EXIT_FAILURE);
    }
  }
  return layers;
}

MultilayerMetadata parse_multilayer_metadata(std::fstream &file) {
  int line_idx = 0;
  bool has_list_prefix;
  int indent = -1;
  std::string field_name;
  ParsedValue value;
  MultilayerMetadata multilayer = {};
  while (parse_line(file, /*min_indent=*/0, /*is_list=*/false, &indent,
                    &has_list_prefix, &line_idx, &field_name, &value)) {
    // Check if string starts with field name.
    if ((field_name == "use_case")) {
      multilayer.use_case = (MultilayerUseCase)value.IntValueInRange(
          /*min=*/0, /*max=*/63, line_idx);
    } else if ((field_name == "layers")) {
      multilayer.layers =
          parse_multilayer_layer_metadata(file,
                                          /*min_indent=*/indent + 1, &line_idx);
    } else {
      fprintf(stderr, "Error: Unknown field %s at line %d\n",
              field_name.c_str(), line_idx);
      exit(EXIT_FAILURE);
    }
  }
  return multilayer;
}

std::string format_depth_representation_element(
    const std::pair<DepthRepresentationElement, bool> &element) {
  if (!element.second) {
    return "absent";
  } else {
    return std::to_string(
               depth_representation_element_to_double(element.first)) +
           " (sign " + std::to_string(element.first.sign_flag) + " exponent " +
           std::to_string(element.first.exponent) + " mantissa " +
           std::to_string(element.first.mantissa) + " mantissa_len " +
           std::to_string(element.first.mantissa_len) + ")";
    ;
  }
}

std::string format_color_properties(
    const std::pair<ColorProperties, bool> &color_properties) {
  if (!color_properties.second) {
    return "absent";
  } else {
    return std::to_string(color_properties.first.color_primaries) + "/" +
           std::to_string(color_properties.first.transfer_characteristics) +
           "/" + std::to_string(color_properties.first.matrix_coefficients) +
           (color_properties.first.color_range ? "F" : "L");
  }
}

void validate_multilayer_metadata(const MultilayerMetadata &multilayer) {
  if (multilayer.layers.empty()) {
    fprintf(stderr, "Error: No layers found, there must be at least one\n");
    exit(EXIT_FAILURE);
  }
  if (multilayer.layers.size() > 4) {
    fprintf(stderr, "Error: Too many layers, found %d, max 4\n",
            (int)multilayer.layers.size());
    exit(EXIT_FAILURE);
  }

  bool same_view_type = true;
  MultilayerViewType view_type = multilayer.layers[0].layer_view_type;
  for (const LayerMetadata &layer : multilayer.layers) {
    if (layer.layer_view_type != view_type) {
      same_view_type = false;
      break;
    }
  }

  for (int i = 0; i < (int)multilayer.layers.size(); ++i) {
    const LayerMetadata &layer = multilayer.layers[i];
    switch (multilayer.use_case) {
      case MULTILAYER_USE_CASE_GLOBAL_ALPHA:
      case MULTILAYER_USE_CASE_GLOBAL_DEPTH:
      case MULTILAYER_USE_CASE_STEREO:
      case MULTILAYER_USE_CASE_STEREO_GLOBAL_ALPHA:
      case MULTILAYER_USE_CASE_STEREO_GLOBAL_DEPTH:
      case MULTILAYER_USE_CASE_444_GLOBAL_ALPHA:
      case MULTILAYER_USE_CASE_444_GLOBAL_DEPTH:
        if (layer.layer_metadata_scope != SCOPE_GLOBAL) {
          fprintf(
              stderr,
              "Error: for use_case %d, all layers must have scope %d, found %d "
              "instead for layer %d (zero-based index)\n",
              multilayer.use_case, SCOPE_GLOBAL, layer.layer_metadata_scope, i);
          exit(EXIT_FAILURE);
        }
        break;
      default: break;
    }
    switch (multilayer.use_case) {
      case MULTILAYER_USE_CASE_GLOBAL_ALPHA:
      case MULTILAYER_USE_CASE_GLOBAL_DEPTH:
      case MULTILAYER_USE_CASE_ALPHA:
      case MULTILAYER_USE_CASE_DEPTH:
      case MULTILAYER_USE_CASE_444_GLOBAL_ALPHA:
      case MULTILAYER_USE_CASE_444_GLOBAL_DEPTH:
      case MULTILAYER_USE_CASE_444:
      case MULTILAYER_USE_CASE_420_444:
        if (!same_view_type) {
          fprintf(stderr,
                  "Error: for use_case %d, all layers must have the same view "
                  "type, found different view_type for layer %d (zero-based "
                  "index)\n",
                  multilayer.use_case, i);
          exit(EXIT_FAILURE);
        }
      default: break;
    }
    if (layer.layer_type != MULTILAYER_LAYER_TYPE_UNSPECIFIED)
      switch (multilayer.use_case) {
        case MULTILAYER_USE_CASE_GLOBAL_ALPHA:
        case MULTILAYER_USE_CASE_ALPHA:
        case MULTILAYER_USE_CASE_STEREO_GLOBAL_ALPHA:
        case MULTILAYER_USE_CASE_STEREO_ALPHA:
          if (layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE &&
              layer.layer_type != MULTILAYER_LAYER_TYPE_ALPHA) {
            fprintf(stderr,
                    "Error: for use_case %d, all layers must be of type %d or "
                    "%d, found %d for layer %d (zero-based index)\n",
                    multilayer.use_case, MULTILAYER_LAYER_TYPE_TEXTURE,
                    MULTILAYER_LAYER_TYPE_ALPHA, layer.layer_type, i);
            exit(EXIT_FAILURE);
          }
          break;
        case MULTILAYER_USE_CASE_GLOBAL_DEPTH:
        case MULTILAYER_USE_CASE_DEPTH:
        case MULTILAYER_USE_CASE_STEREO_GLOBAL_DEPTH:
        case MULTILAYER_USE_CASE_STEREO_DEPTH:
          if (layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE &&
              layer.layer_type != MULTILAYER_LAYER_TYPE_DEPTH) {
            fprintf(stderr,
                    "Error: for use_case %d, all layers must be of type %d or "
                    "%d, found %d for layer %d (zero-based index)\n",
                    multilayer.use_case, MULTILAYER_LAYER_TYPE_TEXTURE,
                    MULTILAYER_LAYER_TYPE_DEPTH, layer.layer_type, i);
            exit(EXIT_FAILURE);
          }
          break;
        case MULTILAYER_USE_CASE_STEREO:
          if (layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE) {
            fprintf(stderr,
                    "Error: for use_case %d, all layers must be of type %d, "
                    "found %d for layer %d (zero-based index)\n",
                    multilayer.use_case, MULTILAYER_LAYER_TYPE_TEXTURE,
                    layer.layer_type, i);
            exit(EXIT_FAILURE);
          }
          break;
        case MULTILAYER_USE_CASE_444_GLOBAL_ALPHA:
          if (layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE_1 &&
              layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE_2 &&
              layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE_3 &&
              layer.layer_type != MULTILAYER_LAYER_TYPE_ALPHA) {
            fprintf(stderr,
                    "Error: for use_case %d, all layers must be of type %d, "
                    "%d, %d, or %d, found %d for layer %d (zero-based index)\n",
                    multilayer.use_case, MULTILAYER_LAYER_TYPE_TEXTURE_1,
                    MULTILAYER_LAYER_TYPE_TEXTURE_2,
                    MULTILAYER_LAYER_TYPE_TEXTURE_3,
                    MULTILAYER_LAYER_TYPE_ALPHA, layer.layer_type, i);
            exit(EXIT_FAILURE);
          }
          break;
        case MULTILAYER_USE_CASE_444_GLOBAL_DEPTH:
          if (layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE_1 &&
              layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE_2 &&
              layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE_3 &&
              layer.layer_type != MULTILAYER_LAYER_TYPE_DEPTH) {
            fprintf(stderr,
                    "Error: for use_case %d, all layers must be of type %d, "
                    "%d, %d, or %d, found %d for layer %d (zero-based index)\n",
                    multilayer.use_case, MULTILAYER_LAYER_TYPE_TEXTURE_1,
                    MULTILAYER_LAYER_TYPE_TEXTURE_2,
                    MULTILAYER_LAYER_TYPE_TEXTURE_3,
                    MULTILAYER_LAYER_TYPE_DEPTH, layer.layer_type, i);
            exit(EXIT_FAILURE);
          }
          break;
        case MULTILAYER_USE_CASE_444:
          if (layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE_1 &&
              layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE_2 &&
              layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE_3) {
            fprintf(
                stderr,
                "Error: for use_case %d, all layers must be of type %d, %d, or "
                "%d, found %d for layer %d (zero-based index)\n",
                multilayer.use_case, MULTILAYER_LAYER_TYPE_TEXTURE_1,
                MULTILAYER_LAYER_TYPE_TEXTURE_2,
                MULTILAYER_LAYER_TYPE_TEXTURE_3, layer.layer_type, i);
            exit(EXIT_FAILURE);
          }
          break;
        case MULTILAYER_USE_CASE_420_444:
          if (layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE &&
              layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE_1 &&
              layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE_2 &&
              layer.layer_type != MULTILAYER_LAYER_TYPE_TEXTURE_3) {
            fprintf(stderr,
                    "Error: for use_case %d, all layers must be of type %d, "
                    "%d, %d, or %d, found %d for layer %d (zero-based index)\n",
                    multilayer.use_case, MULTILAYER_LAYER_TYPE_TEXTURE,
                    MULTILAYER_LAYER_TYPE_TEXTURE_1,
                    MULTILAYER_LAYER_TYPE_TEXTURE_2,
                    MULTILAYER_LAYER_TYPE_TEXTURE_3, layer.layer_type, i);
            exit(EXIT_FAILURE);
          }
          break;
        default: break;
      }
    if (layer.layer_dependency_idc >= (1 << i)) {
      fprintf(stderr,
              "Error: layer_dependency_idc of layer %d (zero-based index) must "
              "be in [0, %d], found %d for layer %d (zero-based index)\n",
              i, (1 << i) - 1, layer.layer_dependency_idc, i);
      exit(EXIT_FAILURE);
    }
    if ((layer.layer_type == MULTILAYER_LAYER_TYPE_ALPHA ||
         layer.layer_type == MULTILAYER_LAYER_TYPE_DEPTH) &&
        layer.layer_color_description.second) {
      fprintf(stderr,
              "Error: alpha or depth layers cannot have "
              "layer_color_description for layer %d (zero-based index)\n",
              i);
      exit(EXIT_FAILURE);
    }
  }
}

}  // namespace

MultilayerMetadata parse_multilayer_file(const char *metadata_path) {
  std::fstream file(metadata_path);
  if (!file.is_open()) {
    fprintf(stderr, "Error: Failed to open %s\n", metadata_path);
    exit(EXIT_FAILURE);
  }

  const MultilayerMetadata multilayer = parse_multilayer_metadata(file);
  validate_multilayer_metadata(multilayer);
  return multilayer;
}

void print_multilayer_metadata(const MultilayerMetadata &multilayer) {
  printf("=== Multilayer metadata ===\n");
  printf("use_case: %d\n", multilayer.use_case);
  for (size_t i = 0; i < multilayer.layers.size(); ++i) {
    const LayerMetadata &layer = multilayer.layers[i];
    printf("layer %d\n", (int)i);
    printf("  layer_type: %d\n", layer.layer_type);
    printf("  luma_plane_only_flag: %d\n", layer.luma_plane_only_flag);
    printf("  layer_view_type: %d\n", layer.layer_view_type);
    printf("  group_id: %d\n", layer.group_id);
    printf("  layer_dependency_idc: %d\n", layer.layer_dependency_idc);
    printf("  layer_metadata_scope: %d\n", layer.layer_metadata_scope);
    printf("  layer_color_description: %s\n",
           format_color_properties(layer.layer_color_description).c_str());
    if (layer.layer_type == MULTILAYER_LAYER_TYPE_ALPHA) {
      printf("  alpha:\n");
      printf("    alpha_use_idc: %d\n", layer.global_alpha_info.alpha_use_idc);
      printf("    alpha_bit_depth: %d\n",
             layer.global_alpha_info.alpha_bit_depth);
      printf("    alpha_clip_idc: %d\n",
             layer.global_alpha_info.alpha_clip_idc);
      printf("    alpha_incr_flag: %d\n",
             layer.global_alpha_info.alpha_incr_flag);
      printf("    alpha_transparent_value: %hu\n",
             layer.global_alpha_info.alpha_transparent_value);
      printf("    alpha_opaque_value: %hu\n",
             layer.global_alpha_info.alpha_opaque_value);
      printf("    alpha_color_description: %s\n",
             format_color_properties(
                 layer.global_alpha_info.alpha_color_description)
                 .c_str());
      printf("    label_type_id:");
      for (uint16_t label_type_id : layer.global_alpha_info.label_type_id) {
        printf(" %d", label_type_id);
      }
      printf("\n");
    } else if (layer.layer_type == MULTILAYER_LAYER_TYPE_DEPTH) {
      printf("  depth:\n");
      printf("    z_near: %s\n",
             format_depth_representation_element(layer.global_depth_info.z_near)
                 .c_str());
      printf("    z_far: %s\n",
             format_depth_representation_element(layer.global_depth_info.z_far)
                 .c_str());
      printf("    d_min: %s\n",
             format_depth_representation_element(layer.global_depth_info.d_min)
                 .c_str());
      printf("    d_max: %s\n",
             format_depth_representation_element(layer.global_depth_info.d_max)
                 .c_str());
      printf("    depth_representation_type: %d\n",
             layer.global_depth_info.depth_representation_type);
      printf("    disparity_ref_view_id: %d\n",
             layer.global_depth_info.disparity_ref_view_id);
      printf("    depth_nonlinear_precision: %d\n",
             layer.global_depth_info.depth_nonlinear_precision);
      printf("    depth_nonlinear_representation_model:");
      for (uint32_t depth_nonlinear_representation_model :
           layer.global_depth_info.depth_nonlinear_representation_model) {
        printf(" %d", depth_nonlinear_representation_model);
      }
      printf("\n");
    }
  }
  printf("\n");
}

}  // namespace libaom_examples
