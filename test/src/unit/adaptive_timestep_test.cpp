/**
 * @file adaptive_timestep_test.cpp
 * @brief 单元测试：验证 `adaptive_timestep` 的核心语义、边界条件和回归行为。
 */
#include <adaptive_timestep/adaptive_timestep.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{

struct ScriptedSolve
{
    bool converged{true};
    int iterations{5};
    double wall{0.01};
};

class FakeBackend
{
public:
    explicit FakeBackend(std::vector<ScriptedSolve> script)
        : script_(std::move(script))
    {
    }

    double currentTime() const noexcept { return time_; }
    void setTimeStep(double dt) { dt_ = dt; attemptedDts.push_back(dt); }
    void beginAttempt() { ++beginCount; controlChangedThisAttempt_ = false; }

    MPMC::NonlinearSolveResult solve()
    {
        if (solveIndex_ >= script_.size())
            throw std::runtime_error("Fake solve script exhausted.");
        const auto s = script_[solveIndex_++];
        MPMC::NonlinearSolveResult result;
        result.converged = s.converged;
        result.nonlinearIterations = s.iterations;
        result.reasonCode = s.converged ? 2 : -5;
        result.reason = s.converged ? "converged" : "diverged";
        result.wallTimeSeconds = s.wall;
        result.linearIterations = 2LL * s.iterations;
        return result;
    }

    bool updateWellControls()
    {
        if (controlSwitchesRemaining_ > 0 && !controlChangedThisAttempt_)
        {
            --controlSwitchesRemaining_;
            controlChangedThisAttempt_ = true;
            return true;
        }
        return false;
    }

    void rejectAttempt()
    {
        ++rejectCount;
        controlChangedThisAttempt_ = false;
    }

    void acceptAttempt(double time)
    {
        ++acceptCount;
        time_ = time;
        acceptedTimes.push_back(time);
    }

    void alignCurrentTime(double time) { time_ = time; }

    void setControlSwitches(int count) { controlSwitchesRemaining_ = count; }

    std::vector<double> attemptedDts;
    std::vector<double> acceptedTimes;
    int beginCount{0};
    int rejectCount{0};
    int acceptCount{0};

private:
    std::vector<ScriptedSolve> script_;
    std::size_t solveIndex_{0};
    double time_{0.0};
    double dt_{0.0};
    int controlSwitchesRemaining_{0};
    bool controlChangedThisAttempt_{false};
};

void require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

bool near(double a, double b, double tol = 1.0e-12)
{
    return std::abs(a - b) <= tol;
}

void testNonlinearStagnationDetector()
{
    MPMC::NonlinearStagnationConfig cfg;
    cfg.enabled = true;
    cfg.minimumIterations = 8;
    cfg.stagnantIterations = 5;
    cfg.requiredRelativeImprovement = 1.0e-4;

    MPMC::NonlinearStagnationDetector detector(cfg);

    // Mirrors the failed three-EOS jobs: one useful Newton step, followed by
    // an essentially unchanged residual plateau.  The guard must wait until
    // the configured minimum iteration before declaring divergence.
    const std::vector<double> plateau{
        3.109317,
        1.042624,
        1.058129,
        1.058129,
        1.058129,
        1.058129,
        1.058129,
        1.058129,
        1.058129};
    for (std::size_t i = 0; i + 1 < plateau.size(); ++i)
    {
        require(!detector.update(static_cast<int>(i), plateau[i]),
                "Stagnation guard fired before minimumIterations.");
    }
    require(detector.update(8, plateau.back()),
            "Stagnation guard missed a persistent residual plateau.");

    detector.reset();
    for (int iteration = 0; iteration <= 20; ++iteration)
    {
        const double residual = 3.0 * std::pow(0.85, iteration);
        require(!detector.update(iteration, residual),
                "Stagnation guard terminated an improving nonlinear solve.");
    }

    MPMC::NonlinearStagnationConfig disabledCfg = cfg;
    disabledCfg.enabled = false;
    MPMC::NonlinearStagnationDetector disabled(disabledCfg);
    for (int iteration = 0; iteration <= 20; ++iteration)
    {
        require(!disabled.update(iteration, 1.0),
                "Disabled stagnation guard changed solver behavior.");
    }
}

void testFailureRetryAndExactOutput()
{
    MPMC::AdaptiveTimeStepConfig cfg;
    cfg.fixedOutputDt = 1.0;
    cfg.minimumDt = 0.125;
    cfg.cutFactor = 0.5;
    cfg.growthFactor = 2.0;
    cfg.easyNonlinearIterations = 4;
    cfg.difficultNonlinearIterations = 8;

    FakeBackend backend({
        {false, 8, 0.10},
        {true, 3, 0.05},
        {true, 3, 0.05}});

    MPMC::AdaptiveTimeStepper<FakeBackend> stepper(backend, cfg);

    std::vector<double> outputTimes;
    stepper.run(1, [&](std::size_t, double time) { outputTimes.push_back(time); });

    require(backend.attemptedDts.size() == 3, "Unexpected attempt count.");
    require(near(backend.attemptedDts[0], 1.0), "First dt changed.");
    require(near(backend.attemptedDts[1], 0.5), "Rejected dt did not halve.");
    require(near(backend.attemptedDts[2], 0.5), "Clipped remainder is wrong.");
    require(backend.rejectCount == 1, "Reject count mismatch.");
    require(backend.acceptCount == 2, "Accept count mismatch.");
    require(outputTimes.size() == 1 && near(outputTimes[0], 1.0), "Output time changed.");
    require(near(backend.currentTime(), 1.0), "Final time is not exact.");

    const auto &stats = stepper.statistics();
    require(stats.rejectedSteps == 1, "Statistics rejectedSteps mismatch.");
    require(stats.acceptedSteps == 2, "Statistics acceptedSteps mismatch.");
    require(stats.nonlinearSolves == 3, "Statistics solve count mismatch.");
    require(stats.attemptedLinearIterations == 28, "Attempted KSP statistics mismatch.");
    require(stats.acceptedLinearIterations == 12, "Accepted KSP statistics mismatch.");
    require(near(stats.minimumAcceptedDt, 0.5), "Minimum accepted dt mismatch.");
    require(near(stats.maximumAcceptedDt, 0.5), "Maximum accepted dt mismatch.");
    require(near(stats.totalSolveTimeSeconds, 0.20), "Solve-time accumulation mismatch.");
    require(near(stats.acceptedSolveTimeSeconds, 0.10), "Accepted solve-time mismatch.");
}

void testClippedStepPreservesRecommendation()
{
    MPMC::AdaptiveTimeStepConfig cfg;
    cfg.fixedOutputDt = 1.0;
    cfg.minimumDt = 0.01;
    cfg.growthFactor = 2.0;
    cfg.easyNonlinearIterations = 10;
    cfg.difficultNonlinearIterations = 20;

    FakeBackend backend({
        {false, 30, 0.0},
        {true, 5, 0.0},
        {true, 5, 0.0},
        {true, 5, 0.0},
        {true, 5, 0.0}});

    MPMC::AdaptiveTimeStepper<FakeBackend> stepper(backend, cfg);
    stepper.run(2, [](std::size_t, double) {});

    // Failure: 1 -> 0.5. Easy 0.5 step grows recommendation back to 1.0.
    // The 0.5 remainder to the first output point is clipped from a 1.0
    // recommendation and must not lower it; the next interval therefore starts 1.0.
    require(backend.attemptedDts.size() >= 4, "Not enough attempts for clipped test.");
    require(near(backend.attemptedDts[3], 1.0), "Clipped output step polluted next recommendation.");
}

void testMaximumInternalStepIsIndependentOfOutputCadence()
{
    MPMC::AdaptiveTimeStepConfig cfg;
    cfg.fixedOutputDt = 60.0;
    cfg.maximumDt = 2.0;
    cfg.minimumDt = 0.01;
    cfg.growthFactor = 2.0;
    cfg.easyNonlinearIterations = 10;
    cfg.difficultNonlinearIterations = 20;

    std::vector<ScriptedSolve> script(
        30, ScriptedSolve{true, 3, 0.0});
    FakeBackend backend(std::move(script));
    MPMC::AdaptiveTimeStepper<FakeBackend> stepper(backend, cfg);

    stepper.run(1, [](std::size_t, double) {});

    require(backend.attemptedDts.size() == 30,
            "60 s output interval did not split into 2 s internal steps.");
    for (double dt : backend.attemptedDts)
        require(near(dt, 2.0), "Internal dt exceeded the configured maximum.");
    require(near(backend.currentTime(), 60.0),
            "Internal dt cap changed the mandatory output time.");
}

void testDifficultAcceptedStepShrinks()
{
    MPMC::AdaptiveTimeStepConfig cfg;
    cfg.fixedOutputDt = 2.0;
    cfg.minimumDt = 0.1;
    cfg.easyNonlinearIterations = 4;
    cfg.difficultNonlinearIterations = 8;
    cfg.difficultShrinkFactor = 0.8;

    FakeBackend backend({
        {true, 9, 0.0},
        {true, 5, 0.0},
        {true, 5, 0.0}});

    MPMC::AdaptiveTimeStepper<FakeBackend> stepper(backend, cfg);
    stepper.run(2, [](std::size_t, double) {});

    require(backend.attemptedDts.size() >= 2, "Difficult-step test missing second attempt.");
    require(near(backend.attemptedDts[1], 1.6), "Difficult accepted step did not shrink by 0.8.");
}

void testControlSwitchResolveCountsAllIterations()
{
    MPMC::AdaptiveTimeStepConfig cfg;
    cfg.fixedOutputDt = 1.0;
    cfg.minimumDt = 0.1;
    cfg.easyNonlinearIterations = 4;
    cfg.difficultNonlinearIterations = 8;
    cfg.maximumWellControlIterations = 3;

    FakeBackend backend({
        {true, 3, 0.0},
        {true, 6, 0.0}});
    backend.setControlSwitches(1);

    MPMC::AdaptiveTimeStepper<FakeBackend> stepper(backend, cfg);
    stepper.run(1, [](std::size_t, double) {});

    const auto &stats = stepper.statistics();
    require(stats.controlSwitches == 1, "Control switch not counted.");
    require(stats.nonlinearSolves == 2, "Control re-solve not counted.");
    require(stats.acceptedNonlinearIterations == 9, "Accepted step did not accumulate control-cycle iterations.");
}

void testMinimumDtFailure()
{
    MPMC::AdaptiveTimeStepConfig cfg;
    cfg.fixedOutputDt = 1.0;
    cfg.minimumDt = 0.6;
    cfg.cutFactor = 0.5;

    FakeBackend backend({{false, 1, 0.0}});
    MPMC::AdaptiveTimeStepper<FakeBackend> stepper(backend, cfg);

    bool threw = false;
    try
    {
        stepper.run(1, [](std::size_t, double) {});
    }
    catch (const MPMC::AdaptiveTimeStepFailure &)
    {
        threw = true;
    }

    require(threw, "minimumDt failure was not reported.");
    require(backend.rejectCount == 1, "Rejected state was not restored before failure.");
}

} // namespace

int main()
{
    testNonlinearStagnationDetector();
    testFailureRetryAndExactOutput();
    testClippedStepPreservesRecommendation();
    testMaximumInternalStepIsIndependentOfOutputCadence();
    testDifficultAcceptedStepShrinks();
    testControlSwitchResolveCountsAllIterations();
    testMinimumDtFailure();

    std::cout << "AdaptiveTimeStepper core validation: ALL PASS\n";
    return 0;
}
