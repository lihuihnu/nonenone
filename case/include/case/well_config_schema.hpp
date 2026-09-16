/**
 * @file well_config_schema.hpp
 * @brief 正式算例共享的井配置数据结构；物理参数仍由各 case 的 well_config.hpp 保存。
 */
#pragma once

#include <cstddef>
#include <limits>
#include <vector>

namespace MPMC::cases::well_config
{

enum class Type { Injector, Producer };
enum class Control { Bhp, TotalRate, OilRate, GasRate, WaterRate };
enum class InjectionPhase { Oil, Gas, Water };

template <int DefaultKCount = -1>
struct Completion
{
    int i{0};
    int j{0};
    int kBegin{0};
    int kCount{DefaultKCount};
};

struct Perforation
{
    std::size_t inputCell{0};
    double wellIndex{0.0};
};

/** Config must expose Config::InitialState::pressure. */
template <class Config>
struct StructuredWellDefinition
{
    int id{0};
    const char *name{""};
    Type type{Type::Producer};
    Control control{Control::Bhp};
    double target{0.0};
    double initialBhp{Config::InitialState::pressure};
    Completion<> completion{};
    double radius{0.10};
    double skin{0.0};
    InjectionPhase injectionPhase{InjectionPhase::Gas};
    int injectedComponent{-1};
};

template <class Config, class CompletionType>
struct ScheduledStructuredDefaults
{
    static constexpr Control control = Control::Bhp;
    inline static constexpr CompletionType completion{};
    static constexpr int injectedComponent = -1;
    static constexpr double minimumBhp = std::numeric_limits<double>::quiet_NaN();
    static constexpr double maximumBhp = std::numeric_limits<double>::quiet_NaN();
    static constexpr bool scheduleEnabled = true;
    static constexpr double openTime = 0.0;
    static constexpr double closeTime = std::numeric_limits<double>::max();
};

/** Field order preserves the former scheduled structured-case aggregate API. */
template <class Config,
          class CompletionType = Completion<>,
          class Defaults = ScheduledStructuredDefaults<Config, CompletionType>>
struct ScheduledStructuredWellDefinition
{
    int id{0};
    const char *name{""};
    Type type{Type::Producer};
    Control control{Defaults::control};
    double target{0.0};
    double initialBhp{Config::InitialState::pressure};
    CompletionType completion{Defaults::completion};
    double radius{0.10};
    double skin{0.0};
    InjectionPhase injectionPhase{InjectionPhase::Gas};
    int injectedComponent{Defaults::injectedComponent};
    double minimumBhp{Defaults::minimumBhp};
    double maximumBhp{Defaults::maximumBhp};
    bool scheduleEnabled{Defaults::scheduleEnabled};
    double openTime{Defaults::openTime};
    double closeTime{Defaults::closeTime};
};

template <class Config>
struct PerforatedDefaults
{
    static constexpr double radius = 0.10;
};

/** Field order preserves the former explicit-perforation aggregate API. */
template <class Config, class Defaults = PerforatedDefaults<Config>>
struct PerforatedWellDefinition
{
    int id{0};
    const char *name{""};
    Type type{Type::Producer};
    Control control{Control::Bhp};
    double target{0.0};
    double initialBhp{Config::InitialState::pressure};
    std::vector<Perforation> perforations;
    double radius{Defaults::radius};
    double skin{0.0};
    InjectionPhase injectionPhase{InjectionPhase::Gas};
    int injectedComponent{-1};
    double minimumBhp{std::numeric_limits<double>::quiet_NaN()};
    double maximumBhp{std::numeric_limits<double>::quiet_NaN()};
    double maximumWaterRate{std::numeric_limits<double>::quiet_NaN()};
    bool scheduleEnabled{true};
    double openTime{0.0};
    double closeTime{std::numeric_limits<double>::max()};
};

} // namespace MPMC::cases::well_config
