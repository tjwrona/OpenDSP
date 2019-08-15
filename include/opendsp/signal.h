#pragma once

#include <complex>
#include <vector>

namespace opendsp
{

template <typename T>
class Signal final
{
public:
    /** @brief Support the use of iterators */
    using iterator = typename std::vector<T>::iterator;

    /** @brief Support the use of const_iterators */
    using const_iterator = typename std::vector<T>::const_iterator;

    Signal(const std::size_t sampleRate, const std::size_t size)
        : sampleRate(sampleRate)
        , samples(size)
    {}

    Signal(const std::size_t sampleRate, const std::vector<T>& samples)
        : sampleRate(sampleRate)
        , samples(samples)
    {}

    Signal(const Signal& other)
        : sampleRate(other.sampleRate)
        , samples(other.samples)
    {}

    Signal(Signal&& other)
        : sampleRate(std::move(other.sampleRate))
        , samples(std::move(other.samples))
    {}

    Signal& operator=(Signal other) noexcept
    {
        this->swap(other);
        return *this;
    }

    Signal& operator=(Signal&& other)
    {
        this->swap(other);
        return *this;
    }

    ~Signal() = default;

    std::size_t GetSampleRate() const
    {
        return sampleRate;
    }

    std::size_t GetLength() const
    {
        return samples.size();
    }

    iterator begin()
    {
        return samples.begin();
    }

    const_iterator begin() const
    {
        return samples.begin();
    }

    const_iterator cbegin() const
    {
        return samples.cbegin();
    }

    iterator end()
    {
        return samples.end();
    }

    const_iterator end() const
    {
        return samples.end();
    }

    const_iterator cend() const
    {
        return samples.cend();
    }

    T& operator[](const std::size_t index)
    {
        return samples[index];
    }

    const T& operator[](const std::size_t index) const
    {
        return samples[index];
    }

    void swap(Signal& other)
    {
        using std::swap;
        swap(sampleRate, other.sampleRate);
        swap(samples, other.samples);
    }

private:
    std::size_t sampleRate;
    std::vector<T> samples;
};

template <typename T>
void swap(Signal<T>& a, Signal<T>& b)
{
    a.swap(b);
}

} /* namespace opendsp */