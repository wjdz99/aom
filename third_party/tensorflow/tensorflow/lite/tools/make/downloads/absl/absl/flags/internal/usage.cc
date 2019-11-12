//
// Copyright 2019 The Abseil Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "absl/flags/internal/usage.h"

#include <map>
#include <string>

#include "absl/flags/flag.h"
#include "absl/flags/internal/path_util.h"
#include "absl/flags/internal/program_name.h"
#include "absl/flags/usage_config.h"
#include "absl/strings/ascii.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_split.h"
#include "absl/synchronization/mutex.h"

ABSL_FLAG(bool, help, false,
          "show help on important flags for this binary [tip: all flags can "
          "have two dashes]");
ABSL_FLAG(bool, helpfull, false, "show help on all flags");
ABSL_FLAG(bool, helpshort, false,
          "show help on only the main module for this program");
ABSL_FLAG(bool, helppackage, false,
          "show help on all modules in the main package");
ABSL_FLAG(bool, version, false, "show version and build info and exit");
ABSL_FLAG(bool, only_check_args, false, "exit after checking all flags");
ABSL_FLAG(std::string, helpon, "",
          "show help on the modules named by this flag value");
ABSL_FLAG(std::string, helpmatch, "",
          "show help on modules whose name contains the specified substr");

namespace absl {
namespace flags_internal {
namespace {

// This class is used to emit an XML element with `tag` and `text`.
// It adds opening and closing tags and escapes special characters in the text.
// For example:
// std::cout << XMLElement("title", "Milk & Cookies");
// prints "<title>Milk &amp; Cookies</title>"
class XMLElement {
 public:
  XMLElement(absl::string_view tag, absl::string_view txt)
      : tag_(tag), txt_(txt) {}

  friend std::ostream& operator<<(std::ostream& out,
                                  const XMLElement& xml_elem) {
    out << "<" << xml_elem.tag_ << ">";

    for (auto c : xml_elem.txt_) {
      switch (c) {
        case '"':
          out << "&quot;";
          break;
        case '\'':
          out << "&apos;";
          break;
        case '&':
          out << "&amp;";
          break;
        case '<':
          out << "&lt;";
          break;
        case '>':
          out << "&gt;";
          break;
        default:
          out << c;
          break;
      }
    }

    return out << "</" << xml_elem.tag_ << ">";
  }

 private:
  absl::string_view tag_;
  absl::string_view txt_;
};

// --------------------------------------------------------------------
// Helper class to pretty-print info about a flag.

class FlagHelpPrettyPrinter {
 public:
  // Pretty printer holds on to the std::ostream& reference to direct an output
  // to that stream.
  FlagHelpPrettyPrinter(int max_line_len, std::ostream* out)
      : out_(*out),
        max_line_len_(max_line_len),
        line_len_(0),
        first_line_(true) {}

  void Write(absl::string_view str, bool wrap_line = false) {
    // Empty std::string - do nothing.
    if (str.empty()) return;

    std::vector<absl::string_view> tokens;
    if (wrap_line) {
      tokens = absl::StrSplit(str, absl::ByAnyChar(" \f\n\r\t\v"),
                              absl::SkipEmpty());
    } else {
      tokens.push_back(str);
    }

    for (auto token : tokens) {
      bool new_line = (line_len_ == 0);

      // Write the token, ending the std::string first if necessary/possible.
      if (!new_line && (line_len_ + token.size() >= max_line_len_)) {
        EndLine();
        new_line = true;
      }

      if (new_line) {
        StartLine();
      } else {
        out_ << ' ';
        ++line_len_;
      }

      out_ << token;
      line_len_ += token.size();
    }
  }

  void StartLine() {
    if (first_line_) {
      out_ << "    ";
      line_len_ = 4;
      first_line_ = false;
    } else {
      out_ << "      ";
      line_len_ = 6;
    }
  }
  void EndLine() {
    out_ << '\n';
    line_len_ = 0;
  }

 private:
  std::ostream& out_;
  const int max_line_len_;
  int line_len_;
  bool first_line_;
};

void FlagHelpHumanReadable(const flags_internal::CommandLineFlag& flag,
                           std::ostream* out) {
  FlagHelpPrettyPrinter printer(80, out);  // Max line length is 80.

  // Flag name.
  printer.Write(absl::StrCat("-", flag.Name()));

  // Flag help.
  printer.Write(absl::StrCat("(", flag.Help(), ");"), /*wrap_line=*/true);

  // Flag data type (for V1 flags only).
  if (!flag.IsAbseilFlag() && !flag.IsRetired()) {
    printer.Write(absl::StrCat("type: ", flag.Typename(), ";"));
  }

  // The listed default value will be the actual default from the flag
  // definition in the originating source file, unless the value has
  // subsequently been modified using SetCommandLineOption() with mode
  // SET_FLAGS_DEFAULT.
  std::string dflt_val = flag.DefaultValue();
  if (flag.IsOfType<std::string>()) {
    dflt_val = absl::StrCat("\"", dflt_val, "\"");
  }
  printer.Write(absl::StrCat("default: ", dflt_val, ";"));

  if (flag.modified) {
    std::string curr_val = flag.CurrentValue();
    if (flag.IsOfType<std::string>()) {
      curr_val = absl::StrCat("\"", curr_val, "\"");
    }
    printer.Write(absl::StrCat("currently: ", curr_val, ";"));
  }

  printer.EndLine();
}

// Shows help for every filename which matches any of the filters
// If filters are empty, shows help for every file.
// If a flag's help message has been stripped (e.g. by adding '#define
// STRIP_FLAG_HELP 1' then this flag will not be displayed by '--help'
// and its variants.
void FlagsHelpImpl(std::ostream& out, flags_internal::FlagKindFilter filter_cb,
                   HelpFormat format = HelpFormat::kHumanReadable) {
  if (format == HelpFormat::kHumanReadable) {
    out << flags_internal::ShortProgramInvocationName() << ": "
        << flags_internal::ProgramUsageMessage() << "\n\n";
  } else {
    // XML schema is not a part of our public API for now.
    out << "<?xml version=\"1.0\"?>\n"
        // The document.
        << "<AllFlags>\n"
        // The program name and usage.
        << XMLElement("program", flags_internal::ShortProgramInvocationName())
        << '\n'
        << XMLElement("usage", flags_internal::ProgramUsageMessage()) << '\n';
  }

  // Map of package name to
  //   map of file name to
  //     vector of flags in the file.
  // This map is used to output matching flags grouped by package and file
  // name.
  std::map<std::string,
           std::map<std::string,
                    std::vector<const flags_internal::CommandLineFlag*>>>
      matching_flags;

  flags_internal::ForEachFlag([&](flags_internal::CommandLineFlag* flag) {
    absl::MutexLock l(InitFlagIfNecessary(flag));

    std::string flag_filename = flag->Filename();

    // Ignore retired flags.
    if (flag->IsRetired()) return;

    // If the flag has been stripped, pretend that it doesn't exist.
    if (flag->Help() == flags_internal::kStrippedFlagHelp) return;

    // Make sure flag satisfies the filter
    if (!filter_cb || !filter_cb(flag_filename)) return;

    matching_flags[std::string(flags_internal::Package(flag_filename))]
                  [flag_filename]
                      .push_back(flag);
  });

  absl::string_view
      package_separator;             // controls blank lines between packages.
  absl::string_view file_separator;  // controls blank lines between files.
  for (const auto& package : matching_flags) {
    if (format == HelpFormat::kHumanReadable) {
      out << package_separator;
      package_separator = "\n\n";
    }

    file_separator = "";
    for (const auto& flags_in_file : package.second) {
      if (format == HelpFormat::kHumanReadable) {
        out << file_separator << "  Flags from " << flags_in_file.first
            << ":\n";
        file_separator = "\n";
      }

      for (const auto* flag : flags_in_file.second) {
        flags_internal::FlagHelp(out, *flag, format);
      }
    }
  }

  if (format == HelpFormat::kHumanReadable) {
    if (filter_cb && matching_flags.empty()) {
      out << "  No modules matched: use -helpfull\n";
    }
  } else {
    // The end of the document.
    out << "</AllFlags>\n";
  }
}

ABSL_CONST_INIT absl::Mutex usage_message_guard(absl::kConstInit);
ABSL_CONST_INIT std::string* program_usage_message
    GUARDED_BY(usage_message_guard) = nullptr;

}  // namespace

// --------------------------------------------------------------------
// Sets the "usage" message to be used by help reporting routines.

void SetProgramUsageMessage(absl::string_view new_usage_message) {
  absl::MutexLock l(&usage_message_guard);

  if (flags_internal::program_usage_message != nullptr) {
    ABSL_INTERNAL_LOG(FATAL, "SetProgramUsageMessage() called twice.");
    std::exit(1);
  }

  program_usage_message = new std::string(new_usage_message);
}

// --------------------------------------------------------------------
// Returns the usage message set by SetProgramUsageMessage().
// Note: We able to return string_view here only because calling
// SetProgramUsageMessage twice is prohibited.
absl::string_view ProgramUsageMessage() {
  absl::MutexLock l(&usage_message_guard);

  return program_usage_message != nullptr
             ? absl::string_view(*program_usage_message)
             : "Warning: SetProgramUsageMessage() never called";
}

// --------------------------------------------------------------------
// Produces the help message describing specific flag.
void FlagHelp(std::ostream& out, const flags_internal::CommandLineFlag& flag,
              HelpFormat format) {
  if (format == HelpFormat::kHumanReadable)
    flags_internal::FlagHelpHumanReadable(flag, &out);
}

// --------------------------------------------------------------------
// Produces the help messages for all flags matching the filter.
// If filter is empty produces help messages for all flags.
void FlagsHelp(std::ostream& out, absl::string_view filter, HelpFormat format) {
  flags_internal::FlagKindFilter filter_cb = [&](absl::string_view filename) {
    return filter.empty() || filename.find(filter) != absl::string_view::npos;
  };
  flags_internal::FlagsHelpImpl(out, filter_cb, format);
}

// --------------------------------------------------------------------
// Checks all the 'usage' command line flags to see if any have been set.
// If so, handles them appropriately.
int HandleUsageFlags(std::ostream& out) {
  if (absl::GetFlag(FLAGS_helpshort)) {
    flags_internal::FlagsHelpImpl(
        out, flags_internal::GetUsageConfig().contains_helpshort_flags,
        HelpFormat::kHumanReadable);
    return 1;
  }

  if (absl::GetFlag(FLAGS_helpfull)) {
    // show all options
    flags_internal::FlagsHelp(out);
    return 1;
  }

  if (!absl::GetFlag(FLAGS_helpon).empty()) {
    flags_internal::FlagsHelp(
        out, absl::StrCat("/", absl::GetFlag(FLAGS_helpon), "."));
    return 1;
  }

  if (!absl::GetFlag(FLAGS_helpmatch).empty()) {
    flags_internal::FlagsHelp(out, absl::GetFlag(FLAGS_helpmatch));
    return 1;
  }

  if (absl::GetFlag(FLAGS_help)) {
    flags_internal::FlagsHelpImpl(
        out, flags_internal::GetUsageConfig().contains_help_flags);

    out << "\nTry --helpfull to get a list of all flags.\n";

    return 1;
  }

  if (absl::GetFlag(FLAGS_helppackage)) {
    flags_internal::FlagsHelpImpl(
        out, flags_internal::GetUsageConfig().contains_helppackage_flags);

    out << "\nTry --helpfull to get a list of all flags.\n";

    return 1;
  }

  if (absl::GetFlag(FLAGS_version)) {
    if (flags_internal::GetUsageConfig().version_string)
      out << flags_internal::GetUsageConfig().version_string();
    // Unlike help, we may be asking for version in a script, so return 0
    return 0;
  }

  if (absl::GetFlag(FLAGS_only_check_args)) {
    return 0;
  }

  return -1;
}

}  // namespace flags_internal
}  // namespace absl
