use anyhow::{Context, Result};

/// A single hourly weather record.
#[derive(Debug, Clone)]
pub struct HourlyRecord {
    /// Month (1-12).
    pub month: u8,
    /// Day of month (1-31).
    pub day: u8,
    /// Hour (1-24).
    pub hour: u8,
    /// Dry bulb temperature in °C.
    pub dry_bulb_temperature: f64,
    /// Relative humidity in %.
    pub relative_humidity: f64,
    /// Global horizontal radiation in Wh/m^2.
    pub global_horizontal_radiation: f64,
    /// Direct normal radiation in Wh/m^2.
    pub direct_normal_radiation: f64,
    /// Diffuse horizontal radiation in Wh/m^2.
    pub diffuse_horizontal_radiation: f64,
    /// Wind speed in m/s.
    pub wind_speed: f64,
    /// Wind direction in degrees from north.
    pub wind_direction: f64,
}

/// Parsed EPW weather data.
#[derive(Debug, Clone)]
pub struct WeatherData {
    /// Location name.
    pub location: String,
    /// Latitude in degrees.
    pub latitude: f64,
    /// Longitude in degrees.
    pub longitude: f64,
    /// Time zone (hours from UTC).
    pub timezone: f64,
    /// Elevation in meters.
    pub elevation: f64,
    /// 8760 hourly records (or 8784 for leap years).
    pub records: Vec<HourlyRecord>,
}

impl WeatherData {
    /// Parses EPW (EnergyPlus Weather) file content.
    ///
    /// EPW format: 8 header lines followed by hourly data rows.
    /// Each data row has 35 fields, comma-separated.
    pub fn from_epw(content: &str) -> Result<Self> {
        let lines: Vec<&str> = content.lines().collect();
        if lines.len() < 9 {
            anyhow::bail!("EPW file too short: expected at least 9 lines");
        }

        // Parse LOCATION header (line 0)
        // Format: LOCATION,city,state_province,country,source,WMO,lat,lon,tz,elevation
        let location_fields: Vec<&str> = lines[0].split(',').collect();
        if location_fields.len() < 10 {
            anyhow::bail!("Invalid LOCATION header");
        }

        let location = format!(
            "{}, {}",
            location_fields[1].trim(),
            location_fields[3].trim()
        );
        let latitude: f64 = location_fields[6]
            .trim()
            .parse()
            .context("Invalid latitude")?;
        let longitude: f64 = location_fields[7]
            .trim()
            .parse()
            .context("Invalid longitude")?;
        let timezone: f64 = location_fields[8]
            .trim()
            .parse()
            .context("Invalid timezone")?;
        let elevation: f64 = location_fields[9]
            .trim()
            .parse()
            .context("Invalid elevation")?;

        // Data rows start at line 8
        let mut records = Vec::new();
        for (i, line) in lines.iter().enumerate().skip(8) {
            let fields: Vec<&str> = line.split(',').collect();
            if fields.len() < 35 {
                continue; // Skip malformed lines
            }

            let record = HourlyRecord {
                month: fields[1]
                    .trim()
                    .parse()
                    .with_context(|| format!("Invalid month at line {i}"))?,
                day: fields[2]
                    .trim()
                    .parse()
                    .with_context(|| format!("Invalid day at line {i}"))?,
                hour: fields[3]
                    .trim()
                    .parse()
                    .with_context(|| format!("Invalid hour at line {i}"))?,
                dry_bulb_temperature: fields[6]
                    .trim()
                    .parse()
                    .with_context(|| format!("Invalid dry bulb at line {i}"))?,
                relative_humidity: fields[8]
                    .trim()
                    .parse()
                    .with_context(|| format!("Invalid RH at line {i}"))?,
                global_horizontal_radiation: fields[13]
                    .trim()
                    .parse()
                    .with_context(|| format!("Invalid GHR at line {i}"))?,
                direct_normal_radiation: fields[14]
                    .trim()
                    .parse()
                    .with_context(|| format!("Invalid DNR at line {i}"))?,
                diffuse_horizontal_radiation: fields[15]
                    .trim()
                    .parse()
                    .with_context(|| format!("Invalid DHR at line {i}"))?,
                wind_speed: fields[21]
                    .trim()
                    .parse()
                    .with_context(|| format!("Invalid wind speed at line {i}"))?,
                wind_direction: fields[20]
                    .trim()
                    .parse()
                    .with_context(|| format!("Invalid wind dir at line {i}"))?,
            };
            records.push(record);
        }

        Ok(Self {
            location,
            latitude,
            longitude,
            timezone,
            elevation,
            records,
        })
    }

    /// Creates simple synthetic weather data for testing.
    ///
    /// Generates 8760 hours with sinusoidal temperature variation.
    pub fn synthetic(
        location: &str,
        latitude: f64,
        longitude: f64,
        mean_temp: f64,
        temp_amplitude: f64,
    ) -> Self {
        let mut records = Vec::with_capacity(8760);
        let days_in_month: [u16; 12] = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

        for month_idx in 0..12 {
            for day in 1..=days_in_month[month_idx] {
                for hour in 1..=24u16 {
                    let day_of_year: f64 =
                        days_in_month[..month_idx].iter().sum::<u16>() as f64 + day as f64;
                    let hour_of_year = (day_of_year - 1.0) * 24.0 + hour as f64;

                    // Annual sinusoidal temperature variation (peak in July/August)
                    let annual_phase = 2.0 * std::f64::consts::PI * (day_of_year - 200.0) / 365.0;
                    // Daily variation (peak at 14:00)
                    let daily_phase = 2.0 * std::f64::consts::PI * (hour as f64 - 14.0) / 24.0;

                    let temp =
                        mean_temp + temp_amplitude * annual_phase.cos() + 3.0 * daily_phase.cos();

                    // Simple solar radiation model
                    let solar_altitude_factor = if (7..=19).contains(&hour) {
                        let solar_hour = (hour as f64 - 12.0) / 6.0;
                        (1.0 - solar_hour * solar_hour).max(0.0)
                    } else {
                        0.0
                    };
                    let ghr = 800.0 * solar_altitude_factor;

                    records.push(HourlyRecord {
                        month: (month_idx + 1) as u8,
                        day: day as u8,
                        hour: hour as u8,
                        dry_bulb_temperature: temp,
                        relative_humidity: 60.0,
                        global_horizontal_radiation: ghr,
                        direct_normal_radiation: ghr * 0.6,
                        diffuse_horizontal_radiation: ghr * 0.4,
                        wind_speed: 3.0,
                        wind_direction: 180.0,
                    });

                    let _ = hour_of_year; // suppress unused
                }
            }
        }

        Self {
            location: location.to_string(),
            latitude,
            longitude,
            timezone: 0.0,
            elevation: 0.0,
            records,
        }
    }

    /// Returns the number of hours in the dataset.
    pub fn num_hours(&self) -> usize {
        self.records.len()
    }

    /// Returns the annual mean temperature.
    pub fn mean_temperature(&self) -> f64 {
        if self.records.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.records.iter().map(|r| r.dry_bulb_temperature).sum();
        sum / self.records.len() as f64
    }

    /// Returns the design heating temperature (99th percentile cold).
    pub fn design_heating_temperature(&self) -> f64 {
        if self.records.is_empty() {
            return 0.0;
        }
        let mut temps: Vec<f64> = self
            .records
            .iter()
            .map(|r| r.dry_bulb_temperature)
            .collect();
        temps.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let idx = (temps.len() as f64 * 0.01) as usize;
        temps[idx]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_synthetic_weather() {
        let weather = WeatherData::synthetic("Test City", 52.0, 13.0, 10.0, 12.0);
        assert_eq!(weather.num_hours(), 8760);
        assert_eq!(weather.location, "Test City");

        let mean = weather.mean_temperature();
        assert!(
            (mean - 10.0).abs() < 3.0,
            "Mean temp should be near 10°C, got {mean}"
        );
    }

    #[test]
    fn test_design_heating_temperature() {
        let weather = WeatherData::synthetic("Test City", 52.0, 13.0, 10.0, 12.0);
        let design_temp = weather.design_heating_temperature();
        assert!(
            design_temp < 0.0,
            "Design heating temp should be below 0°C, got {design_temp}"
        );
    }

    #[test]
    fn test_epw_parse_minimal() {
        // Minimal EPW-like content
        let header = "LOCATION,TestCity,State,Country,Source,123456,52.0,13.0,1.0,50.0\n\
                       DESIGN CONDITIONS,0\n\
                       TYPICAL/EXTREME PERIODS,0\n\
                       GROUND TEMPERATURES,0\n\
                       HOLIDAYS/DAYLIGHT SAVINGS,No,0,0,0\n\
                       COMMENTS 1,test\n\
                       COMMENTS 2,test\n\
                       DATA PERIODS,1,1,Data,Sunday,1/1,12/31\n";

        // Build one data line with 35 fields
        // Fields: year,month,day,hour,minute,source,drybulb,dewpoint,relhum,atmpressure,
        //         exthoriz,extdirect,horizinfra,ghr,dnr,dhr,...(remaining filled with 0)
        let data_line = "2020,1,1,1,60,?,5.0,2.0,80,101325,0,0,0,0,0,0,0,0,0,0,180,3.0,0,0,0,0,0,0,0,0,0,0,0,0,0\n";

        let content = format!("{}{}", header, data_line);
        let weather = WeatherData::from_epw(&content).unwrap();

        assert_eq!(weather.location, "TestCity, Country");
        assert!((weather.latitude - 52.0).abs() < 1e-10);
        assert_eq!(weather.records.len(), 1);
        assert!((weather.records[0].dry_bulb_temperature - 5.0).abs() < 1e-10);
    }
}
