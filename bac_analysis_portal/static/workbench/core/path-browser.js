export function createPathBrowser(openBrowser) {
  const targets = new Map();

  return {
    open(selector, options = {}) {
      return openBrowser(selector, options);
    },
    registerTarget(selector, handler) {
      if (!selector || typeof handler !== "function") return;
      targets.set(selector, handler);
    },
    async handleSelection(selector, context) {
      const handler = targets.get(selector);
      if (!handler) return false;
      await handler(context);
      return true;
    },
  };
}
