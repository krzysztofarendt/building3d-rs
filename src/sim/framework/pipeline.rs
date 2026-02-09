use anyhow::Result;

use super::{Bus, SimContext, SimModule};

/// Executes a sequence of simulation modules.
pub struct Pipeline {
    modules: Vec<Box<dyn SimModule>>,
}

impl Pipeline {
    pub fn new() -> Self {
        Self { modules: vec![] }
    }

    pub fn with_module<M: SimModule + 'static>(mut self, module: M) -> Self {
        self.modules.push(Box::new(module));
        self
    }

    pub fn init(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        for module in self.modules.iter_mut() {
            module.init(ctx, bus)?;
        }
        Ok(())
    }

    pub fn step(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        for module in self.modules.iter_mut() {
            module.step(ctx, bus)?;
        }
        Ok(())
    }
}

impl Default for Pipeline {
    fn default() -> Self {
        Self::new()
    }
}
