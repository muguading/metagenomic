export function createNavigation() {
  const modules = new Map();

  return {
    register(tabId, domainModule) {
      if (!tabId || !domainModule) return;
      modules.set(tabId, domainModule);
    },
    async activate(tabId) {
      await modules.get(tabId)?.activate?.();
    },
    async refresh(tabId) {
      await modules.get(tabId)?.refresh?.();
    },
  };
}
