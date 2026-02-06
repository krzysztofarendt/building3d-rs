/// A time-varying schedule for internal gains, setpoints, etc.
///
/// Provides hourly values that repeat on a weekly or daily basis.
#[derive(Debug, Clone)]
pub struct Schedule {
    pub name: String,
    /// Hourly values for the schedule period.
    /// If 24 values: repeats daily.
    /// If 168 values: repeats weekly (Mon-Sun, 24h each).
    /// If 8760 values: annual (no repeat).
    values: Vec<f64>,
}

impl Schedule {
    /// Creates a schedule from a list of hourly values.
    pub fn new(name: &str, values: Vec<f64>) -> Self {
        Self {
            name: name.to_string(),
            values,
        }
    }

    /// Creates a constant schedule.
    pub fn constant(name: &str, value: f64) -> Self {
        Self {
            name: name.to_string(),
            values: vec![value],
        }
    }

    /// Creates a typical office occupancy schedule (8am-6pm weekdays).
    pub fn office_occupancy() -> Self {
        let mut values = Vec::with_capacity(168);
        for day in 0..7 {
            for hour in 0..24 {
                let is_weekday = day < 5;
                let is_working_hour = (8..18).contains(&hour);
                let val = if is_weekday && is_working_hour {
                    1.0
                } else {
                    0.0
                };
                values.push(val);
            }
        }
        Self::new("office_occupancy", values)
    }

    /// Creates a residential occupancy schedule.
    pub fn residential_occupancy() -> Self {
        let mut values = Vec::with_capacity(24);
        for hour in 0..24 {
            let val = match hour {
                0..=6 => 1.0,   // sleeping
                7..=8 => 0.5,   // morning
                9..=16 => 0.2,  // away
                17..=21 => 0.8, // evening
                22..=23 => 1.0, // night
                _ => 0.0,
            };
            values.push(val);
        }
        Self::new("residential_occupancy", values)
    }

    /// Gets the schedule value for a given hour of the year (0-8759).
    pub fn value_at(&self, hour_of_year: usize) -> f64 {
        if self.values.is_empty() {
            return 0.0;
        }
        if self.values.len() == 1 {
            return self.values[0];
        }
        let idx = hour_of_year % self.values.len();
        self.values[idx]
    }
}

/// Internal gains profile combining occupancy, equipment, and lighting.
#[derive(Debug, Clone)]
pub struct InternalGainsProfile {
    /// Heat gain per person in W.
    pub heat_per_person: f64,
    /// Number of occupants.
    pub num_occupants: f64,
    /// Occupancy schedule.
    pub occupancy: Schedule,
    /// Equipment power in W.
    pub equipment_power: f64,
    /// Equipment schedule.
    pub equipment: Schedule,
    /// Lighting power in W.
    pub lighting_power: f64,
    /// Lighting schedule.
    pub lighting: Schedule,
}

impl InternalGainsProfile {
    /// Creates a default office gains profile.
    pub fn office(floor_area: f64) -> Self {
        Self {
            heat_per_person: 120.0, // seated, light office work
            num_occupants: floor_area / 10.0,
            occupancy: Schedule::office_occupancy(),
            equipment_power: floor_area * 15.0, // 15 W/m^2
            equipment: Schedule::office_occupancy(),
            lighting_power: floor_area * 10.0, // 10 W/m^2
            lighting: Schedule::office_occupancy(),
        }
    }

    /// Total internal gains at a given hour of the year.
    pub fn gains_at(&self, hour_of_year: usize) -> f64 {
        let occ = self.occupancy.value_at(hour_of_year);
        let eq = self.equipment.value_at(hour_of_year);
        let lit = self.lighting.value_at(hour_of_year);

        self.heat_per_person * self.num_occupants * occ
            + self.equipment_power * eq
            + self.lighting_power * lit
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constant_schedule() {
        let s = Schedule::constant("test", 42.0);
        assert!((s.value_at(0) - 42.0).abs() < 1e-10);
        assert!((s.value_at(5000) - 42.0).abs() < 1e-10);
    }

    #[test]
    fn test_office_occupancy() {
        let s = Schedule::office_occupancy();
        assert_eq!(s.values.len(), 168);

        // Monday 10am (hour index = 0*24+10 = 10)
        assert!(
            (s.value_at(10) - 1.0).abs() < 1e-10,
            "Office occupied Mon 10am"
        );

        // Monday 2am (hour index = 0*24+2 = 2)
        assert!((s.value_at(2) - 0.0).abs() < 1e-10, "Office empty Mon 2am");

        // Saturday 10am (hour index = 5*24+10 = 130)
        assert!(
            (s.value_at(130) - 0.0).abs() < 1e-10,
            "Office empty Saturday"
        );
    }

    #[test]
    fn test_residential_occupancy() {
        let s = Schedule::residential_occupancy();
        assert_eq!(s.values.len(), 24);

        // 3am
        assert!((s.value_at(3) - 1.0).abs() < 1e-10, "Home occupied at 3am");
        // 12pm
        assert!((s.value_at(12) - 0.2).abs() < 1e-10, "Mostly away at noon");
        // 8pm
        assert!((s.value_at(20) - 0.8).abs() < 1e-10, "Home in evening");
    }

    #[test]
    fn test_internal_gains_profile() {
        let profile = InternalGainsProfile::office(100.0); // 100 m^2
        // Monday 10am (hour 10 of year)
        let gains = profile.gains_at(10);
        // Occupancy=1.0: 120*10*1.0 + 1500*1.0 + 1000*1.0 = 1200+1500+1000 = 3700 W
        assert!(
            (gains - 3700.0).abs() < 1.0,
            "Expected ~3700W office gains at 10am, got {gains}"
        );
    }
}
