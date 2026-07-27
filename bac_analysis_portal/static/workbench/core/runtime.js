export function createWorkbenchRuntime({ api, ui, navigation, pathBrowser, constants = {}, selectors = {}, actions = {} }) {
  return Object.freeze({ api, ui, navigation, pathBrowser, constants, selectors, actions });
}
