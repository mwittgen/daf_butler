# This file is part of daf_butler.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import annotations

"""Support for Storage Classes."""

__all__ = ("StorageClass", "StorageClassFactory", "StorageClassConfig")

import builtins
import copy
import itertools
import logging
from typing import Any, Collection, Dict, List, Mapping, Optional, Sequence, Set, Tuple, Type, Union

from lsst.utils import doImportType
from lsst.utils.classes import Singleton
from lsst.utils.introspection import get_full_type_name

from .config import Config, ConfigSubset
from .configSupport import LookupKey
from .storageClassDelegate import StorageClassDelegate

log = logging.getLogger(__name__)


class StorageClassConfig(ConfigSubset):
    """Configuration class for defining Storage Classes."""

    component = "storageClasses"
    defaultConfigFile = "storageClasses.yaml"


class StorageClass:
    """Class describing how a label maps to a particular Python type.

    Parameters
    ----------
    name : `str`
        Name to use for this class.
    pytype : `type` or `str`
        Python type (or name of type) to associate with the `StorageClass`
    components : `dict`, optional
        `dict` mapping name of a component to another `StorageClass`.
    derivedComponents : `dict`, optional
        `dict` mapping name of a derived component to another `StorageClass`.
    parameters : `~collections.abc.Sequence` or `~collections.abc.Set`
        Parameters understood by this `StorageClass` that can control
        reading of data from datastores.
    delegate : `str`, optional
        Fully qualified name of class supporting assembly and disassembly
        of a `pytype` instance.
    converters : `dict` [`str`, `str`], optional
        Mapping of python type to function that can be called to convert
        that python type to the valid type of this storage class.
    """

    _cls_name: str = "BaseStorageClass"
    _cls_components: Optional[Dict[str, StorageClass]] = None
    _cls_derivedComponents: Optional[Dict[str, StorageClass]] = None
    _cls_parameters: Optional[Union[Set[str], Sequence[str]]] = None
    _cls_delegate: Optional[str] = None
    _cls_pytype: Optional[Union[Type, str]] = None
    _cls_converters: Optional[Dict[str, str]] = None
    defaultDelegate: Type = StorageClassDelegate
    defaultDelegateName: str = get_full_type_name(defaultDelegate)

    def __init__(
        self,
        name: Optional[str] = None,
        pytype: Optional[Union[Type, str]] = None,
        components: Optional[Dict[str, StorageClass]] = None,
        derivedComponents: Optional[Dict[str, StorageClass]] = None,
        parameters: Optional[Union[Sequence, Set]] = None,
        delegate: Optional[str] = None,
        converters: Optional[Dict[str, str]] = None,
    ):
        if name is None:
            name = self._cls_name
        if pytype is None:
            pytype = self._cls_pytype
        if components is None:
            components = self._cls_components
        if derivedComponents is None:
            derivedComponents = self._cls_derivedComponents
        if parameters is None:
            parameters = self._cls_parameters
        if delegate is None:
            delegate = self._cls_delegate

        # Merge converters with class defaults.
        self._converters = {}
        if self._cls_converters is not None:
            self._converters.update(self._cls_converters)
        if converters:
            self._converters.update(converters)

        # Version of converters where the python types have been
        # Do not try to import anything until needed.
        self._converters_by_type: Optional[Dict[Type, Type]] = None

        self.name = name

        if pytype is None:
            pytype = object

        self._pytype: Optional[Type]
        if not isinstance(pytype, str):
            # Already have a type so store it and get the name
            self._pytypeName = get_full_type_name(pytype)
            self._pytype = pytype
        else:
            # Store the type name and defer loading of type
            self._pytypeName = pytype
            self._pytype = None

        if components is not None:
            if len(components) == 1:
                raise ValueError(
                    f"Composite storage class {name} is not allowed to have"
                    f" only one component '{next(iter(components))}'."
                    " Did you mean it to be a derived component?"
                )
            self._components = components
        else:
            self._components = {}
        self._derivedComponents = derivedComponents if derivedComponents is not None else {}
        self._parameters = frozenset(parameters) if parameters is not None else frozenset()
        # if the delegate is not None also set it and clear the default
        # delegate
        self._delegate: Optional[Type]
        self._delegateClassName: Optional[str]
        if delegate is not None:
            self._delegateClassName = delegate
            self._delegate = None
        elif components is not None:
            # We set a default delegate for composites so that a class is
            # guaranteed to support something if it is a composite.
            log.debug("Setting default delegate for %s", self.name)
            self._delegate = self.defaultDelegate
            self._delegateClassName = self.defaultDelegateName
        else:
            self._delegate = None
            self._delegateClassName = None

    @property
    def components(self) -> Dict[str, StorageClass]:
        """Return the components associated with this `StorageClass`."""
        return self._components

    @property
    def derivedComponents(self) -> Dict[str, StorageClass]:
        """Return derived components associated with `StorageClass`."""
        return self._derivedComponents

    @property
    def converters(self) -> Dict[str, str]:
        """Return the type converters supported by this `StorageClass`."""
        return self._converters

    @property
    def converters_by_type(self) -> Dict[Type, Type]:
        """Return the type converters as python types."""
        if self._converters_by_type is None:
            self._converters_by_type = {}

            # Loop over list because the dict can be edited in loop.
            for candidate_type_str, converter_str in list(self.converters.items()):
                if hasattr(builtins, candidate_type_str):
                    candidate_type = getattr(builtins, candidate_type_str)
                else:
                    try:
                        candidate_type = doImportType(candidate_type_str)
                    except ImportError as e:
                        log.warning(
                            "Unable to import type %s associated with storage class %s (%s)",
                            candidate_type_str,
                            self.name,
                            e,
                        )
                        del self.converters[candidate_type_str]
                        continue

                try:
                    converter = doImportType(converter_str)
                except ImportError as e:
                    log.warning(
                        "Unable to import conversion function %s associated with storage class %s "
                        "required to convert type %s (%s)",
                        converter_str,
                        self.name,
                        candidate_type_str,
                        e,
                    )
                    del self.converters[candidate_type_str]
                    continue
                if not callable(converter):
                    # doImportType is annotated to return a Type but in actual
                    # fact it can return Any except ModuleType because package
                    # variables can be accessed. This make mypy believe it
                    # is impossible for the return value to not be a callable
                    # so we must ignore the warning.
                    log.warning(  # type: ignore
                        "Conversion function %s associated with storage class "
                        "%s to convert type %s is not a callable.",
                        converter_str,
                        self.name,
                        candidate_type_str,
                    )
                    del self.converters[candidate_type_str]
                    continue
                self._converters_by_type[candidate_type] = converter
        return self._converters_by_type

    @property
    def parameters(self) -> Set[str]:
        """Return `set` of names of supported parameters."""
        return set(self._parameters)

    @property
    def pytype(self) -> Type:
        """Return Python type associated with this `StorageClass`."""
        if self._pytype is not None:
            return self._pytype

        if hasattr(builtins, self._pytypeName):
            pytype = getattr(builtins, self._pytypeName)
        else:
            pytype = doImportType(self._pytypeName)
        self._pytype = pytype
        return self._pytype

    @property
    def delegateClass(self) -> Optional[Type]:
        """Class to use to delegate type-specific actions."""
        if self._delegate is not None:
            return self._delegate
        if self._delegateClassName is None:
            return None
        delegate_class = doImportType(self._delegateClassName)
        self._delegate = delegate_class
        return self._delegate

    def allComponents(self) -> Mapping[str, StorageClass]:
        """Return all defined components.

        This mapping includes all the derived and read/write components
        for the corresponding storage class.

        Returns
        -------
        comp : `dict` of [`str`, `StorageClass`]
            The component name to storage class mapping.
        """
        components = copy.copy(self.components)
        components.update(self.derivedComponents)
        return components

    def delegate(self) -> StorageClassDelegate:
        """Return an instance of a storage class delegate.

        Returns
        -------
        delegate : `StorageClassDelegate`
            Instance of the delegate associated with this `StorageClass`.
            The delegate is constructed with this `StorageClass`.

        Raises
        ------
        TypeError
            This StorageClass has no associated delegate.
        """
        cls = self.delegateClass
        if cls is None:
            raise TypeError(f"No delegate class is associated with StorageClass {self.name}")
        return cls(storageClass=self)

    def isComposite(self) -> bool:
        """Return Boolean indicating whether this is a composite or not.

        Returns
        -------
        isComposite : `bool`
            `True` if this `StorageClass` is a composite, `False`
            otherwise.
        """
        if self.components:
            return True
        return False

    def _lookupNames(self) -> Tuple[LookupKey, ...]:
        """Keys to use when looking up this DatasetRef in a configuration.

        The names are returned in order of priority.

        Returns
        -------
        names : `tuple` of `LookupKey`
            Tuple of a `LookupKey` using the `StorageClass` name.
        """
        return (LookupKey(name=self.name),)

    def knownParameters(self) -> Set[str]:
        """Return set of all parameters known to this `StorageClass`.

        The set includes parameters understood by components of a composite.

        Returns
        -------
        known : `set`
            All parameter keys of this `StorageClass` and the component
            storage classes.
        """
        known = set(self._parameters)
        for sc in self.components.values():
            known.update(sc.knownParameters())
        return known

    def validateParameters(self, parameters: Collection = None) -> None:
        """Check that the parameters are known to this `StorageClass`.

        Does not check the values.

        Parameters
        ----------
        parameters : `~collections.abc.Collection`, optional
            Collection containing the parameters. Can be `dict`-like or
            `set`-like.  The parameter values are not checked.
            If no parameters are supplied, always returns without error.

        Raises
        ------
        KeyError
            Some parameters are not understood by this `StorageClass`.
        """
        # No parameters is always okay
        if not parameters:
            return

        # Extract the important information into a set. Works for dict and
        # list.
        external = set(parameters)

        diff = external - self.knownParameters()
        if diff:
            s = "s" if len(diff) > 1 else ""
            unknown = "', '".join(diff)
            raise KeyError(f"Parameter{s} '{unknown}' not understood by StorageClass {self.name}")

    def filterParameters(self, parameters: Mapping[str, Any], subset: Collection = None) -> Mapping[str, Any]:
        """Filter out parameters that are not known to this `StorageClass`.

        Parameters
        ----------
        parameters : `Mapping`, optional
            Candidate parameters. Can be `None` if no parameters have
            been provided.
        subset : `~collections.abc.Collection`, optional
            Subset of supported parameters that the caller is interested
            in using.  The subset must be known to the `StorageClass`
            if specified. If `None` the supplied parameters will all
            be checked, else only the keys in this set will be checked.

        Returns
        -------
        filtered : `Mapping`
            Valid parameters. Empty `dict` if none are suitable.

        Raises
        ------
        ValueError
            Raised if the provided subset is not a subset of the supported
            parameters or if it is an empty set.
        """
        if not parameters:
            return {}

        known = self.knownParameters()

        if subset is not None:
            if not subset:
                raise ValueError("Specified a parameter subset but it was empty")
            subset = set(subset)
            if not subset.issubset(known):
                raise ValueError(f"Requested subset ({subset}) is not a subset of known parameters ({known})")
            wanted = subset
        else:
            wanted = known

        return {k: parameters[k] for k in wanted if k in parameters}

    def validateInstance(self, instance: Any) -> bool:
        """Check that the supplied Python object has the expected Python type.

        Parameters
        ----------
        instance : `object`
            Object to check.

        Returns
        -------
        isOk : `bool`
            True if the supplied instance object can be handled by this
            `StorageClass`, False otherwise.
        """
        return isinstance(instance, self.pytype)

    def can_convert(self, other: StorageClass) -> bool:
        """Return `True` if this storage class can convert python types
        in the other storage class.

        Parameters
        ----------
        other : `StorageClass`
            The storage class to check.

        Returns
        -------
        can : `bool`
            `True` if this storage class has a registered converter for
            the python type associated with the other storage class. That
            converter will convert the other python type to the one associated
            with this storage class.
        """
        if other.name == self.name:
            # Identical storage classes are compatible.
            return True

        # It may be that the storage class being compared is not
        # available because the python type can't be imported. In that
        # case conversion must be impossible.
        try:
            other_pytype = other.pytype
        except Exception:
            return False

        # Or even this storage class itself can not have the type imported.
        try:
            self_pytype = self.pytype
        except Exception:
            return False

        if issubclass(other_pytype, self_pytype):
            # Storage classes have different names but the same python type.
            return True

        for candidate_type in self.converters_by_type:
            if issubclass(other_pytype, candidate_type):
                return True
        return False

    def coerce_type(self, incorrect: Any) -> Any:
        """Coerce the supplied incorrect instance to the python type
        associated with this `StorageClass`.

        Parameters
        ----------
        incorrect : `object`
            An object that might be the incorrect type.

        Returns
        -------
        correct : `object`
            An object that matches the python type of this `StorageClass`.
            Can be the same object as given. If `None`, `None` will be
            returned.

        Raises
        ------
        TypeError
            Raised if no conversion can be found.
        """
        if incorrect is None:
            return None

        # Possible this is the correct type already.
        if self.validateInstance(incorrect):
            return incorrect

        # Check each registered converter.
        for candidate_type, converter in self.converters_by_type.items():
            if isinstance(incorrect, candidate_type):
                try:
                    return converter(incorrect)
                except Exception:
                    log.error(
                        "Converter %s failed to convert type %s",
                        get_full_type_name(converter),
                        get_full_type_name(incorrect),
                    )
                    raise
        raise TypeError(
            "Type does not match and no valid converter found to convert"
            f" '{type(incorrect)}' to '{self.pytype}'"
        )

    def __eq__(self, other: Any) -> bool:
        """Equality checks name, pytype name, delegate name, and components."""
        if not isinstance(other, StorageClass):
            return NotImplemented

        if self.name != other.name:
            return False

        # We must compare pytype and delegate by name since we do not want
        # to trigger an import of external module code here
        if self._delegateClassName != other._delegateClassName:
            return False
        if self._pytypeName != other._pytypeName:
            return False

        # Ensure we have the same component keys in each
        if set(self.components.keys()) != set(other.components.keys()):
            return False

        # Same parameters
        if self.parameters != other.parameters:
            return False

        # Ensure that all the components have the same type
        for k in self.components:
            if self.components[k] != other.components[k]:
                return False

        # If we got to this point everything checks out
        return True

    def __hash__(self) -> int:
        return hash(self.name)

    def __repr__(self) -> str:
        optionals: Dict[str, Any] = {}
        if self._pytypeName != "object":
            optionals["pytype"] = self._pytypeName
        if self._delegateClassName is not None:
            optionals["delegate"] = self._delegateClassName
        if self._parameters:
            optionals["parameters"] = self._parameters
        if self.components:
            optionals["components"] = self.components
        if self.converters:
            optionals["converters"] = self.converters

        # order is preserved in the dict
        options = ", ".join(f"{k}={v!r}" for k, v in optionals.items())

        # Start with mandatory fields
        r = f"{self.__class__.__name__}({self.name!r}"
        if options:
            r = r + ", " + options
        r = r + ")"
        return r

    def __str__(self) -> str:
        return self.name


class StorageClassFactory(metaclass=Singleton):
    """Factory for `StorageClass` instances.

    This class is a singleton, with each instance sharing the pool of
    StorageClasses. Since code can not know whether it is the first
    time the instance has been created, the constructor takes no arguments.
    To populate the factory with storage classes, a call to
    `~StorageClassFactory.addFromConfig()` should be made.

    Parameters
    ----------
    config : `StorageClassConfig` or `str`, optional
        Load configuration. In a ButlerConfig` the relevant configuration
        is located in the ``storageClasses`` section.
    """

    def __init__(self, config: Optional[Union[StorageClassConfig, str]] = None):
        self._storageClasses: Dict[str, StorageClass] = {}
        self._configs: List[StorageClassConfig] = []

        # Always seed with the default config
        self.addFromConfig(StorageClassConfig())

        if config is not None:
            self.addFromConfig(config)

    def __str__(self) -> str:
        """Return summary of factory.

        Returns
        -------
        summary : `str`
            Summary of the factory status.
        """
        sep = "\n"
        return f"""Number of registered StorageClasses: {len(self._storageClasses)}

StorageClasses
--------------
{sep.join(f"{s}: {self._storageClasses[s]}" for s in self._storageClasses)}
"""

    def __contains__(self, storageClassOrName: Union[StorageClass, str]) -> bool:
        """Indicate whether the storage class exists in the factory.

        Parameters
        ----------
        storageClassOrName : `str` or `StorageClass`
            If `str` is given existence of the named StorageClass
            in the factory is checked. If `StorageClass` is given
            existence and equality are checked.

        Returns
        -------
        in : `bool`
            True if the supplied string is present, or if the supplied
            `StorageClass` is present and identical.

        Notes
        -----
        The two different checks (one for "key" and one for "value") based on
        the type of the given argument mean that it is possible for
        StorageClass.name to be in the factory but StorageClass to not be
        in the factory.
        """
        if isinstance(storageClassOrName, str):
            return storageClassOrName in self._storageClasses
        elif isinstance(storageClassOrName, StorageClass):
            if storageClassOrName.name in self._storageClasses:
                return storageClassOrName == self._storageClasses[storageClassOrName.name]
        return False

    def addFromConfig(self, config: Union[StorageClassConfig, Config, str]) -> None:
        """Add more `StorageClass` definitions from a config file.

        Parameters
        ----------
        config : `StorageClassConfig`, `Config` or `str`
            Storage class configuration. Can contain a ``storageClasses``
            key if part of a global configuration.
        """
        sconfig = StorageClassConfig(config)
        self._configs.append(sconfig)

        # Since we can not assume that we will get definitions of
        # components or parents before their classes are defined
        # we have a helper function that we can call recursively
        # to extract definitions from the configuration.
        def processStorageClass(name: str, sconfig: StorageClassConfig, msg: str = "") -> None:
            # Maybe we've already processed this through recursion
            if name not in sconfig:
                return
            info = sconfig.pop(name)

            # Always create the storage class so we can ensure that
            # we are not trying to overwrite with a different definition
            components = None

            # Extract scalar items from dict that are needed for
            # StorageClass Constructor
            storageClassKwargs = {k: info[k] for k in ("pytype", "delegate", "parameters") if k in info}

            if "converters" in info:
                storageClassKwargs["converters"] = info["converters"].toDict()

            for compName in ("components", "derivedComponents"):
                if compName not in info:
                    continue
                components = {}
                for cname, ctype in info[compName].items():
                    if ctype not in self:
                        processStorageClass(ctype, sconfig, msg)
                    components[cname] = self.getStorageClass(ctype)

                # Fill in other items
                storageClassKwargs[compName] = components

            # Create the new storage class and register it
            baseClass = None
            if "inheritsFrom" in info:
                baseName = info["inheritsFrom"]
                if baseName not in self:
                    processStorageClass(baseName, sconfig, msg)
                baseClass = type(self.getStorageClass(baseName))

            newStorageClassType = self.makeNewStorageClass(name, baseClass, **storageClassKwargs)
            newStorageClass = newStorageClassType()
            self.registerStorageClass(newStorageClass, msg=msg)

        # In case there is a problem, construct a context message for any
        # error reporting.
        files = [str(f) for f in itertools.chain([sconfig.configFile], sconfig.filesRead) if f]
        context = f"when adding definitions from {', '.join(files)}" if files else ""
        log.debug("Adding definitions from config %s", ", ".join(files))

        for name in list(sconfig.keys()):
            processStorageClass(name, sconfig, context)

    @staticmethod
    def makeNewStorageClass(
        name: str, baseClass: Optional[Type[StorageClass]] = StorageClass, **kwargs: Any
    ) -> Type[StorageClass]:
        """Create a new Python class as a subclass of `StorageClass`.

        Parameters
        ----------
        name : `str`
            Name to use for this class.
        baseClass : `type`, optional
            Base class for this `StorageClass`. Must be either `StorageClass`
            or a subclass of `StorageClass`. If `None`, `StorageClass` will
            be used.

        Returns
        -------
        newtype : `type` subclass of `StorageClass`
            Newly created Python type.
        """
        if baseClass is None:
            baseClass = StorageClass
        if not issubclass(baseClass, StorageClass):
            raise ValueError(f"Base class must be a StorageClass not {baseClass}")

        # convert the arguments to use different internal names
        clsargs = {f"_cls_{k}": v for k, v in kwargs.items() if v is not None}
        clsargs["_cls_name"] = name

        # Some container items need to merge with the base class values
        # so that a child can inherit but override one bit.
        # lists (which you get from configs) are treated as sets for this to
        # work consistently.
        for k in ("components", "parameters", "derivedComponents", "converters"):
            classKey = f"_cls_{k}"
            if classKey in clsargs:
                baseValue = getattr(baseClass, classKey, None)
                if baseValue is not None:
                    currentValue = clsargs[classKey]
                    if isinstance(currentValue, dict):
                        newValue = baseValue.copy()
                    else:
                        newValue = set(baseValue)
                    newValue.update(currentValue)
                    clsargs[classKey] = newValue

        # If we have parameters they should be a frozen set so that the
        # parameters in the class can not be modified.
        pk = "_cls_parameters"
        if pk in clsargs:
            clsargs[pk] = frozenset(clsargs[pk])

        return type(f"StorageClass{name}", (baseClass,), clsargs)

    def getStorageClass(self, storageClassName: str) -> StorageClass:
        """Get a StorageClass instance associated with the supplied name.

        Parameters
        ----------
        storageClassName : `str`
            Name of the storage class to retrieve.

        Returns
        -------
        instance : `StorageClass`
            Instance of the correct `StorageClass`.

        Raises
        ------
        KeyError
            The requested storage class name is not registered.
        """
        return self._storageClasses[storageClassName]

    def registerStorageClass(self, storageClass: StorageClass, msg: Optional[str] = None) -> None:
        """Store the `StorageClass` in the factory.

        Will be indexed by `StorageClass.name` and will return instances
        of the supplied `StorageClass`.

        Parameters
        ----------
        storageClass : `StorageClass`
            Type of the Python `StorageClass` to register.
        msg : `str`, optional
            Additional message string to be included in any error message.

        Raises
        ------
        ValueError
            If a storage class has already been registered with
            that storage class name and the previous definition differs.
        """
        if storageClass.name in self._storageClasses:
            existing = self.getStorageClass(storageClass.name)
            if existing != storageClass:
                errmsg = f" {msg}" if msg else ""
                raise ValueError(
                    f"New definition for StorageClass {storageClass.name} ({storageClass!r}) "
                    f"differs from current definition ({existing!r}){errmsg}"
                )
        else:
            self._storageClasses[storageClass.name] = storageClass

    def _unregisterStorageClass(self, storageClassName: str) -> None:
        """Remove the named StorageClass from the factory.

        Parameters
        ----------
        storageClassName : `str`
            Name of storage class to remove.

        Raises
        ------
        KeyError
            The named storage class is not registered.

        Notes
        -----
        This method is intended to simplify testing of StorageClassFactory
        functionality and it is not expected to be required for normal usage.
        """
        del self._storageClasses[storageClassName]
